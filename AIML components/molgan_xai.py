import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics.pairwise import cosine_similarity
import shap
import os
import gc  # Garbage collection
import matplotlib.pyplot as plt

# Detect correct SMILES column name dynamically
def get_smiles_column(csv_file):
    try:
        df = pd.read_csv(csv_file, nrows=1)  # Read just the first row
        for col in df.columns:
            if col.lower().strip() in ["smiles", "smile", "canonical_smiles"]:  # Common variations
                return col
    except Exception as e:
        print(f" Error reading CSV: {e}")
    return None

# Load molecular data efficiently
def load_molecular_data(csv_file, num_samples=100):
    molecule_vectors = []
    smiles_list = []
    count = 0
    
    smiles_column = get_smiles_column(csv_file)
    if not smiles_column:
        print(" No valid 'smiles' column found in CSV. Check file format.")
        return [], []

    try:
        for chunk in pd.read_csv(csv_file, chunksize=1000, usecols=[smiles_column], on_bad_lines='skip'):
            gc.collect()  # Free memory

            for smiles in chunk[smiles_column].dropna().values:  # Process all valid SMILES in chunk
                if count >= num_samples:
                    break

                mol_vector = smiles_to_vector(smiles)
                if mol_vector is not None:
                    molecule_vectors.append(mol_vector)
                    smiles_list.append(smiles)
                    count += 1

                    if count % 10 == 0:
                        print(f" Processed {count}/{num_samples} molecules")

            if count >= num_samples:
                break

            del chunk  # Delete chunk to free memory
            gc.collect()

    except Exception as e:
        print(f" Error during file processing: {e}")
        return [], []

    print(f" Successfully loaded {len(molecule_vectors)} molecular vectors")
    return np.array(molecule_vectors), smiles_list

# Convert SMILES to molecular fingerprint vector
def smiles_to_vector(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=64), dtype=np.float32)
    except Exception as e:
        print(f" Error processing SMILES: {smiles} → {e}")
    return None

# Generator Model
class MolGenerator(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(MolGenerator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim*2),
            nn.ReLU(),
            nn.Linear(hidden_dim*2, output_dim),
            nn.Sigmoid()  # Changed to Sigmoid for fingerprint-like output
        )

    def forward(self, z):
        return self.model(z)

# Discriminator Model
class MolDiscriminator(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(MolDiscriminator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim//2, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.model(x)

# Find the closest real molecule to the generated vector
def vector_to_closest_smiles(generated_vector, dataset_vectors, dataset_smiles):
    """Finds the closest real molecule to the generated vector"""
    similarities = cosine_similarity([generated_vector], dataset_vectors)[0]
    closest_idx = np.argmax(similarities)
    similarity_score = similarities[closest_idx]
    return dataset_smiles[closest_idx], similarity_score

# Wrapper for PyTorch model that returns NumPy arrays
def model_predict_wrapper(discriminator):
    def predict(x):
        with torch.no_grad():
            x_tensor = torch.tensor(x, dtype=torch.float32)
            return discriminator(x_tensor).cpu().numpy()
    return predict

# Function for model interpretation with SHAP
def shap_explain_discriminator(discriminator, generated_vector_np, background_data_np):
    """Create SHAP explanation for a single generated molecule"""
    try:
        # Create a wrapper function that accepts numpy arrays
        predict_fn = model_predict_wrapper(discriminator)
        
        # Create a sample of background data (subset for efficiency)
        background_sample = background_data_np[:50]  
        
        # Initialize the KernelExplainer with numpy arrays
        explainer = shap.KernelExplainer(predict_fn, background_sample)
        
        # Reshape the generated vector to match expected dimensions
        generated_vector_reshaped = generated_vector_np.reshape(1, -1)
        
        # Get SHAP values
        shap_values = explainer.shap_values(generated_vector_reshaped)
        
        return shap_values, explainer
        
    except Exception as e:
        print(f" SHAP error details: {e}")
        import traceback
        traceback.print_exc()
        return None, None

# Function to visualize SHAP values for a single molecule
def visualize_shap_values(shap_values, feature_names=None):
    """Create a simple bar plot of SHAP values for a single molecule"""
    if shap_values is None:
        print(" No SHAP values available to visualize")
        return
    
    # Generate feature names if not provided
    if feature_names is None:
        feature_names = [f"Bit_{i}" for i in range(64)]  # Assuming 64-bit fingerprint
    
    try:
        # Make sure shap_values is a numpy array and properly shaped
        values = np.array(shap_values[0])
        
        # Get absolute importance
        abs_values = np.abs(values)
        
        # Sort indices by importance (highest to lowest)
        sorted_indices = np.argsort(abs_values.flatten())[-10:]  # Top 10
        
        # Get corresponding values and names
        sorted_values = values.flatten()[sorted_indices]
        sorted_names = [feature_names[int(i)] for i in sorted_indices]
        
        # Create visualization
        plt.figure(figsize=(10, 6))
        plt.barh(range(len(sorted_values)), sorted_values, color=['b' if x > 0 else 'r' for x in sorted_values])
        plt.yticks(range(len(sorted_values)), sorted_names)
        plt.xlabel('SHAP Value (impact on discriminator output)')
        plt.title('Top Feature Importance for Generated Molecule')
        plt.tight_layout()
        plt.savefig('molecule_shap_explanation.png')
        print(" SHAP explanation saved as 'molecule_shap_explanation.png'")
        
        # Print top features
        print("\n Top important molecular fingerprint bits:")
        for i, (idx, val) in enumerate(zip(sorted_indices, sorted_values)):
            print(f" {i+1}. {feature_names[int(idx)]}: SHAP value = {val:.4f}")
            
    except Exception as e:
        print(f" Visualization error: {e}")
        import traceback
        traceback.print_exc()

# Simple approach to analyze feature importance directly
def analyze_simple_feature_importance(discriminator, generated_vector):
    """Create a simple feature importance analysis by perturbing each feature"""
    with torch.no_grad():
        # Get the original prediction
        original_pred = discriminator(generated_vector).item()
        
        importance = []
        # Test each feature by setting it to zero
        for i in range(generated_vector.size(1)):
            # Create a copy with one feature zeroed out
            perturbed = generated_vector.clone()
            perturbed[0, i] = 0.0
            
            # Get new prediction
            new_pred = discriminator(perturbed).item()
            
            # Store the difference (importance)
            importance.append(original_pred - new_pred)
        
        return importance

# Main function
def main():
    csv_file = 'SMILES.csv'

    if not os.path.exists(csv_file):
        print(f" Error: File '{csv_file}' not found!")
        return

    print(f" Loading molecular data from {csv_file}...")

    input_dim = 64  # Morgan Fingerprint size
    hidden_dim = 128  # Increased hidden dimension for better representation
    num_samples = 200  # Load more samples for better training
    batch_size = 32
    epochs = 200  # More epochs for better convergence

    try:
        real_molecular_data, smiles_list = load_molecular_data(csv_file, num_samples=num_samples)
        if len(real_molecular_data) == 0:
            print(" No valid molecular data was loaded. Check your CSV file format.")
            return
    except Exception as e:
        print(f" Failed to load data: {e}")
        return

    output_dim = input_dim

    # Initialize the models (Generator and Discriminator for GAN)
    generator = MolGenerator(input_dim, hidden_dim, output_dim)
    discriminator = MolDiscriminator(output_dim, hidden_dim)

    # Set up the loss function and optimizers
    criterion = nn.BCELoss()
    optimizer_G = optim.Adam(generator.parameters(), lr=0.0002, betas=(0.5, 0.999))
    optimizer_D = optim.Adam(discriminator.parameters(), lr=0.0002, betas=(0.5, 0.999))

    # Convert real molecular data to PyTorch tensors
    real_molecular_data_tensor = torch.tensor(real_molecular_data).float()

    print(" Starting training...")
    for epoch in range(epochs):
        gc.collect()  # Free memory

        actual_batch_size = min(batch_size, len(real_molecular_data))
        idx = np.random.randint(0, len(real_molecular_data), actual_batch_size)
        real_molecules = real_molecular_data_tensor[idx]

        # Train Discriminator
        optimizer_D.zero_grad()
        
        real_labels = torch.ones(actual_batch_size, 1)
        fake_labels = torch.zeros(actual_batch_size, 1)
        
        # Real samples
        output_real = discriminator(real_molecules)
        loss_D_real = criterion(output_real, real_labels)
        
        # Fake samples
        noise = torch.randn(actual_batch_size, input_dim)
        fake_molecules = generator(noise)
        output_fake = discriminator(fake_molecules.detach())
        loss_D_fake = criterion(output_fake, fake_labels)
        
        # Total discriminator loss
        loss_D = loss_D_real + loss_D_fake
        loss_D.backward()
        optimizer_D.step()

        # Train Generator
        optimizer_G.zero_grad()
        
        noise = torch.randn(actual_batch_size, input_dim)
        fake_molecules = generator(noise)
        output = discriminator(fake_molecules)
        loss_G = criterion(output, real_labels)
        
        loss_G.backward()
        optimizer_G.step()

        if epoch % 10 == 0:
            print(f" Epoch {epoch}/{epochs}: Loss D = {loss_D.item():.4f}, Loss G = {loss_G.item():.4f}")
            gc.collect()

    print(" Training complete!")
    torch.save(generator.state_dict(), 'generator2.pth')
    torch.save(discriminator.state_dict(), 'discriminator2.pth')

    print("\n Generating a single molecule...")
    with torch.no_grad():
        # Set a fixed seed for reproducibility
        torch.manual_seed(42)
        
        # Generate a single molecule
        noise = torch.randn(1, input_dim)
        generated_vector = generator(noise)
        
        # Convert to numpy for further processing
        generated_vector_np = generated_vector.cpu().numpy()[0]

        # Find the closest match in real molecules
        closest_smiles, similarity = vector_to_closest_smiles(generated_vector_np, real_molecular_data, smiles_list)

        print(f" Generated a single molecule fingerprint")
        print(f" Closest Known Molecule: {closest_smiles}")
        print(f" Similarity Score: {similarity:.4f}")

        # Create molecule object for visualization if needed
        mol = Chem.MolFromSmiles(closest_smiles)
        if mol:
            print(f" Molecule is valid according to RDKit")
        else:
            print(f" Warning: Could not parse SMILES with RDKit")

    # SHAP Explanation for the single generated molecule
    print("\n Generating SHAP explanation for the molecule...")
    try:
        # Convert background data to numpy
        background_data_np = real_molecular_data_tensor[:min(100, len(real_molecular_data_tensor))].cpu().numpy()
        
        # Get SHAP values for the single generated molecule
        shap_values, explainer = shap_explain_discriminator(
            discriminator, 
            generated_vector_np, 
            background_data_np
        )
        
        if shap_values is not None:
            # Create feature names
            feature_names = [f"Bit_{i}" for i in range(input_dim)]
            
            # Visualize and save the SHAP explanation
            visualize_shap_values(shap_values, feature_names)
        else:
            print(" Could not generate SHAP values, using fallback method...")
            
            # Fallback to simple feature importance analysis
            importance = analyze_simple_feature_importance(discriminator, generated_vector)
            
            # Get top features
            feature_names = [f"Bit_{i}" for i in range(input_dim)]
            sorted_indices = np.argsort(np.abs(importance))[-10:]
            
            print("\n Top important molecular fingerprint bits (fallback method):")
            for i, idx in enumerate(reversed(sorted_indices)):
                print(f" {i+1}. {feature_names[idx]}: Importance = {importance[idx]:.4f}")
            
            # Create visualization
            plt.figure(figsize=(10, 6))
            sorted_values = [importance[i] for i in sorted_indices]
            sorted_names = [feature_names[i] for i in sorted_indices]
            plt.barh(range(len(sorted_values)), sorted_values, color=['b' if x > 0 else 'r' for x in sorted_values])
            plt.yticks(range(len(sorted_values)), sorted_names)
            plt.xlabel('Feature Importance (impact on discriminator output)')
            plt.title('Top Feature Importance for Generated Molecule (Fallback Method)')
            plt.tight_layout()
            plt.savefig('molecule_importance_fallback.png')
            print(" Fallback importance explanation saved as 'molecule_importance_fallback.png'")
            
    except Exception as e:
        print(f" Error generating explanations: {e}")
        import traceback
        traceback.print_exc()

    # Clean up
    print("\n Cleaning up resources...")
    del generator, discriminator, real_molecular_data_tensor
    gc.collect()
    print(" Done!")
    
if __name__ == "__main__":
    main()
from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
import torch
import torch.nn as nn
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, QED, rdMolDescriptors, Crippen
import shap
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from typing import Optional, List
import os
import sys
import uvicorn
from pydantic import BaseModel
from contextlib import asynccontextmanager
import logging
import pandas as pd
import pickle

# Add these imports at the top if not already present
import matplotlib
matplotlib.use('Agg')  # Set backend to non-interactive
import matplotlib.pyplot as plt
plt.ioff()  # Turn off interactive mode

# Import the BiomimeticSMILESGenerator class from biome.py
from models.biome import BiomimeticSMILESGenerator

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import the predictor from ml_metrics.py
from models.ml_metrics import DrugDiscoveryMetricsPredictor

# Global variables for models and dataset
generator = None
discriminator = None
real_molecular_data = None
smiles_list = None
biomimetic_model = None

# Model paths
MODEL_DIR = os.path.join(os.path.dirname(__file__), 'models')
GENERATOR_PATH = os.path.join(MODEL_DIR, 'generator2.pth')
DISCRIMINATOR_PATH = os.path.join(MODEL_DIR, 'discriminator2.pth')
BIOMIMETIC_MODEL_PATH = os.path.join(MODEL_DIR, 'biomimetic_model.pkl')

# Initialize the predictor
predictor = DrugDiscoveryMetricsPredictor()

@asynccontextmanager
async def lifespan(app: FastAPI):
    if not initialize_models():
        logger.error("Failed to initialize models")
        raise RuntimeError("Failed to initialize models")
    logger.info("API startup complete")
    yield
    logger.info("Shutting down API")

app = FastAPI(
    title="Molecule Generation API",
    lifespan=lifespan
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Replace "*" with specific origins if needed
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.mount("/static", StaticFiles(directory=MODEL_DIR), name="static")

# Define model classes
class MolGenerator(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int):
        super(MolGenerator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim * 2),
            nn.ReLU(),
            nn.Linear(hidden_dim * 2, output_dim),
            nn.Sigmoid()
        )

    def forward(self, z):
        return self.model(z)

class MolDiscriminator(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int):
        super(MolDiscriminator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim // 2, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.model(x)

def initialize_models():
    global generator, discriminator, biomimetic_model
    try:
        logger.info(f"Current working directory: {os.getcwd()}")
        logger.info(f"Model directory: {MODEL_DIR}")
        
        if not os.path.exists(MODEL_DIR):
            logger.error("Models directory does not exist")
            os.makedirs(MODEL_DIR, exist_ok=True)
            logger.info("Created models directory")
        
        input_dim = 64
        hidden_dim = 128
        generator = MolGenerator(input_dim, hidden_dim, input_dim)
        discriminator = MolDiscriminator(input_dim, hidden_dim)
        
        if os.path.exists(GENERATOR_PATH):
            logger.info(f"Loading generator from {GENERATOR_PATH}")
            generator.load_state_dict(torch.load(GENERATOR_PATH, map_location='cpu'))
            generator.eval()
            logger.info("Generator loaded successfully")
        else:
            logger.error(f"Generator file not found: {GENERATOR_PATH}")
            return False
        
        if os.path.exists(DISCRIMINATOR_PATH):
            logger.info(f"Loading discriminator from {DISCRIMINATOR_PATH}")
            discriminator.load_state_dict(torch.load(DISCRIMINATOR_PATH, map_location='cpu'))
            discriminator.eval()
            logger.info("Discriminator loaded successfully")
        else:
            logger.error(f"Discriminator file not found: {DISCRIMINATOR_PATH}")
            return False

        # Load the biomimetic model using BiomimeticSMILESGenerator
        if os.path.exists(BIOMIMETIC_MODEL_PATH):
            logger.info(f"Loading biomimetic model from {BIOMIMETIC_MODEL_PATH}")
            biomimetic_model = BiomimeticSMILESGenerator()
            biomimetic_model.load_model(BIOMIMETIC_MODEL_PATH)
            logger.info("Biomimetic model loaded successfully")
        else:
            logger.error(f"Biomimetic model file not found: {BIOMIMETIC_MODEL_PATH}")
            raise HTTPException(status_code=500, detail="Biomimetic model not loaded.")
        
        return True
    except Exception as e:
        logger.error(f"Error initializing models: {str(e)}")
        logger.exception("Detailed traceback:")
        return False

@app.get("/api/status")
async def get_status():
    return {
        "status": "online",
        "models_loaded": generator is not None and discriminator is not None,
        "biomimetic_model_loaded": biomimetic_model is not None,
        "data_loaded": real_molecular_data is not None
    }

@app.post("/api/generate")
async def generate_molecule():
    if generator is None or discriminator is None:
        raise HTTPException(status_code=503, detail="Models not initialized")
    
    try:
        with torch.no_grad():
            num_samples = 10
            best_mol = None
            best_score = -1
            best_vector = None
            for _ in range(num_samples):
                noise = torch.randn(1, 64)
                generated_vector = generator(noise)
                generated_vector_np = generated_vector.cpu().numpy()[0]
                disc_score = float(discriminator(generated_vector).item())
                if real_molecular_data is not None and smiles_list is not None:
                    closest_smiles, similarity = vector_to_closest_smiles(
                        generated_vector_np,
                        real_molecular_data,
                        smiles_list
                    )
                    mol = Chem.MolFromSmiles(closest_smiles)
                    if mol and (best_mol is None or disc_score > best_score):
                        best_mol = mol
                        best_score = disc_score
                        best_vector = generated_vector_np
                        best_similarity = similarity
                        best_smiles = closest_smiles

            if best_mol is None:
                raise HTTPException(
                    status_code=500,
                    detail="Failed to generate valid molecule"
                )
            
            result = {
                "vector": best_vector.tolist(),
                "discriminator_score": float(best_score),
                "closest_smiles": best_smiles,
                "similarity_score": float(best_similarity),
                "molecule_image": ""
            }
            img = Draw.MolToImage(best_mol)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            result["molecule_image"] = img_str
            
            logger.info(f"Generated new molecule:")
            logger.info(f"SMILES: {best_smiles}")
            logger.info(f"Discriminator score: {best_score}")
            logger.info(f"Similarity score: {best_similarity}")
            
            return JSONResponse(content=result)
            
    except Exception as e:
        logger.error(f"Error generating molecule: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/generate_biomimetic")
async def generate_biomimetic(target_property_value: float, n_attempts: int = 50):
    """
    Generate a new molecule with the desired property value using the biomimetic model.
    """
    global biomimetic_model
    if biomimetic_model is None:
        raise HTTPException(status_code=500, detail="Biomimetic model not loaded.")
    
    try:
        generated_smiles, predicted_value = biomimetic_model.generate_smiles(
            target_property_value, n_attempts
        )
        if generated_smiles is None:
            raise HTTPException(status_code=400, detail="Failed to generate a molecule.")
        
        return {
            "generated_smiles": generated_smiles,
            "predicted_value": predicted_value,
        }
    except Exception as e:
        logger.error(f"Error generating molecule: {str(e)}")
        raise HTTPException(status_code=500, detail="Error generating molecule.")

def vector_to_closest_smiles(generated_vector, dataset_vectors, dataset_smiles):
    from sklearn.metrics.pairwise import cosine_similarity
    generated_vector = generated_vector.reshape(1, -1)
    dataset_vectors = dataset_vectors.reshape(len(dataset_vectors), -1)
    similarities = cosine_similarity(generated_vector, dataset_vectors)[0]
    top_k = 5
    top_indices = np.argsort(similarities)[-top_k:]
    selected_idx = np.random.choice(top_indices)
    similarity_score = similarities[selected_idx]
    return dataset_smiles[selected_idx], similarity_score

@app.post("/api/upload")
async def upload_file(file: UploadFile = File(...)):
    """
    Upload a dataset and process it for molecule generation.
    """
    global real_molecular_data, smiles_list
    try:
        if not file.filename.endswith('.csv'):
            raise HTTPException(status_code=400, detail="Only CSV files are allowed.")
        
        contents = await file.read()
        df = pd.read_csv(BytesIO(contents))
        
        if 'SMILES' not in df.columns:
            raise HTTPException(status_code=400, detail="CSV must contain a 'SMILES' column.")
        
        smiles_list = df['SMILES'].tolist()
        fingerprints = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=64)
                fingerprints.append(list(fp.ToBitString()))
        
        if not fingerprints:
            raise HTTPException(status_code=400, detail="No valid molecules found in the dataset.")
        
        real_molecular_data = np.array(fingerprints, dtype=np.float32)
        return {"message": "Dataset uploaded successfully.", "molecules_processed": len(smiles_list)}
    except Exception as e:
        logger.error(f"Error uploading file: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Error uploading file: {str(e)}")

@app.post("/api/upload")
async def upload_dataset(file: UploadFile = File(...)):
    """Upload a dataset without processing it."""
    global real_molecular_data, smiles_list
    try:
        if not file.filename.endswith('.csv'):
            raise HTTPException(status_code=400, detail="Only CSV files are allowed.")
        
        contents = await file.read()
        df = pd.read_csv(pd.compat.StringIO(contents.decode('utf-8')))
        
        if 'SMILES' not in df.columns:
            raise HTTPException(status_code=400, detail="CSV must contain a 'SMILES' column.")
        
        smiles_list = df['SMILES'].tolist()
        fingerprints = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=64)
                fingerprints.append(list(fp.ToBitString()))
        
        if not fingerprints:
            raise HTTPException(status_code=400, detail="No valid molecules found in the dataset.")
        
        real_molecular_data = np.array(fingerprints, dtype=np.float32)
        return {"message": "Dataset uploaded successfully.", "molecules_processed": len(smiles_list)}
    except Exception as e:
        logger.error(f"Error uploading dataset: {str(e)}")
        raise HTTPException(status_code=500, detail="Failed to upload dataset.")

class DatasetRequest(BaseModel):
    dataset: str
    model: str

class ProcessRequest(BaseModel):
    model: str  # The selected model (e.g., "normal" or "biomimetic")

@app.post("/api/process")
async def process_dataset(request: ProcessRequest):
    """Process the dataset using the selected model."""
    global real_molecular_data, smiles_list, biomimetic_model, generator, discriminator
    
    logger.info(f"Processing dataset with model: {request.model}")
    
    if real_molecular_data is None or smiles_list is None:
        logger.error("No dataset uploaded.")
        raise HTTPException(status_code=400, detail="No dataset uploaded. Please upload a dataset first.")
    
    model = request.model.strip().lower()
    
    if model == "biomimetic":
        try:
            # Set matplotlib to non-interactive mode
            import matplotlib
            matplotlib.use('Agg')
            plt.ioff()

            target_property_value = 0.8
            max_attempts = 100

            if biomimetic_model is None:
                logger.error("Biomimetic model not loaded.")
                raise HTTPException(status_code=500, detail="Biomimetic model not loaded.")
            
            generated_smiles, predicted_value = None, None
            for _ in range(max_attempts):
                try:
                    # Remove show_plots parameter if generation fails
                    generated_smiles, predicted_value = biomimetic_model.generate_smiles(
                        target_property_value, 
                        n_attempts=1
                    )
                    if generated_smiles is not None:
                        break
                except Exception as e:
                    logger.warning(f"Attempt failed: {str(e)}")
                    continue

            if generated_smiles is None:
                raise HTTPException(status_code=400, detail="Failed to generate valid molecule")

            # Generate molecular structure image
            mol = Chem.MolFromSmiles(generated_smiles)
            if mol is None:
                raise HTTPException(status_code=500, detail="Invalid SMILES generated")

            # Save molecule image
            molecule_image_path = os.path.join(MODEL_DIR, "biomimetic_molecule.png")
            img = Draw.MolToImage(mol)
            img.save(molecule_image_path)
            plt.close('all')  # Close any remaining figures

            response_data = {
                "message": "Dataset processed using the Biomimetic model.",
                "generated_smiles": generated_smiles,
                "similarity_score": float(predicted_value) if predicted_value else 0.0,
                "discriminator_score": 0.0,
                "molecule_image": "/static/biomimetic_molecule.png"
            }
            
            return JSONResponse(content=response_data)

        except Exception as e:
            logger.error(f"Error in biomimetic processing: {str(e)}")
            raise HTTPException(status_code=500, detail=str(e))
            
    elif model == "normal":
        # Use the generator and discriminator models
        if generator is None or discriminator is None:
            logger.error("Generator or discriminator model not loaded.")
            raise HTTPException(status_code=500, detail="Generator or discriminator model not loaded.")
        try:
            with torch.no_grad():
                # Generate molecule
                noise = torch.randn(1, 64)
                generated_vector = generator(noise)
                generated_vector_np = generated_vector.cpu().numpy()[0]
                disc_score = float(discriminator(generated_vector).item())

                # Find closest molecule
                closest_smiles, similarity = vector_to_closest_smiles(
                    generated_vector_np,
                    real_molecular_data,
                    smiles_list
                )
                mol = Chem.MolFromSmiles(closest_smiles)
                if mol is None:
                    raise HTTPException(
                        status_code=500,
                        detail="Failed to generate valid molecule"
                    )

                # Generate SHAP explanation
                class ModelWrapper:
                    def __init__(self, model):
                        self.model = model

                    def __call__(self, x):
                        with torch.no_grad():
                            return self.model(torch.FloatTensor(x)).cpu().numpy()

                wrapped_model = ModelWrapper(discriminator)
                background = real_molecular_data[:100]  # Use first 100 samples as background
                explainer = shap.KernelExplainer(wrapped_model, background)
                shap_values = explainer.shap_values(generated_vector_np.reshape(1, -1))

                # Create SHAP plot
                plt.figure(figsize=(10, 6))
                shap.summary_plot(shap_values, generated_vector_np.reshape(1, -1), show=False)
                shap_graph_path = os.path.join(MODEL_DIR, "shap_explanation.png")
                plt.savefig(shap_graph_path)
                plt.close()

                # Generate molecular structure image
                molecule_image_path = os.path.join(MODEL_DIR, "molecule_image.png")
                img = Draw.MolToImage(mol)
                img.save(molecule_image_path, format="PNG")

                # Convert file paths to URLs
                shap_graph_url = f"/static/{os.path.basename(shap_graph_path)}"
                molecule_image_url = f"/static/{os.path.basename(molecule_image_path)}"

                # Prepare the result
                result = {
                    "message": "Dataset processed using the Normal model.",
                    "generated_smiles": closest_smiles,
                    "discriminator_score": float(disc_score),
                    "similarity_score": float(similarity),
                    "shap_graph": shap_graph_url,
                    "molecule_image": molecule_image_url,
                }
                return JSONResponse(content=result)

        except Exception as e:
            logger.error(f"Error during normal model processing: {str(e)}")
            raise HTTPException(status_code=500, detail=f"Error during normal model processing: {str(e)}")
    
    else:
        logger.error(f"Invalid model selected: {request.model}")
        raise HTTPException(status_code=400, detail="Invalid model selected. Choose 'normal' or 'biomimetic'.")

class MoleculeEvalRequest(BaseModel):
    smiles: str

def format_evaluation_report(results):
    """Format evaluation results into a detailed report"""
    def format_check(value, threshold, higher_is_better=True):
        if higher_is_better:
            return "[PASS]" if value >= threshold else "[FAIL]"
        return "[PASS]" if value <= threshold else "[FAIL]"

    def format_risk(risk):
        risk_symbols = {
            "Low": "[LOW RISK]",
            "Medium": "[MEDIUM RISK]",
            "High": "[HIGH RISK]"
        }
        return f"{risk} {risk_symbols.get(risk, '')}"

    report = f"""
COMPOUND INFORMATION
SMILES: {results['compound_info']['smiles']}
Formula: {results['compound_info']['formula']}
Exact Mass: {results['compound_info']['exact_mass']:.4f} Da
Heavy Atoms: {results['compound_info']['heavy_atoms']}
Ring Count: {results['compound_info']['ring_count']}
Overall QED Score: {results['compound_info']['qed_score']:.3f}

DRUG-LIKENESS (Lipinski's Rule of Five)
Molecular Weight: {results['drug_likeness']['molecular_weight']:.2f} Da {format_check(results['drug_likeness']['molecular_weight'], 500, False)}
LogP: {results['drug_likeness']['logp']:.2f} {format_check(results['drug_likeness']['logp'], 5, False)}
H-Bond Donors: {results['drug_likeness']['hbd']} {format_check(results['drug_likeness']['hbd'], 5, False)}
H-Bond Acceptors: {results['drug_likeness']['hba']} {format_check(results['drug_likeness']['hba'], 10, False)}
Rotatable Bonds: {results['drug_likeness']['rotatable_bonds']} {format_check(results['drug_likeness']['rotatable_bonds'], 10, False)}

Rules Satisfied: {5 - results['drug_likeness']['lipinski_violations']}/5
Overall: {"[DRUG-LIKE]" if results['drug_likeness']['is_drug_like'] else "[NON-DRUG-LIKE]"}

BIOAVAILABILITY
Oral Bioavailability: {results['bioavailability']['oral_bioavailability']:.2f}% {format_check(results['bioavailability']['oral_bioavailability'], 30)}
First-Pass Metabolism: {results['bioavailability']['first_pass_metabolism']:.2f}%
LogD: {results['bioavailability']['logd']:.2f}

Overall: {format_risk(results['bioavailability'].get('overall_rating', 'N/A'))}

TOXICITY PREDICTION
Hepatotoxicity Risk: {format_risk(results['toxicity']['hepatotoxicity'])}
Nephrotoxicity Risk: {format_risk(results['toxicity']['nephrotoxicity'])}
Carcinogenicity Risk: {format_risk(results['toxicity']['carcinogenicity'])}
hERG Inhibition Risk: {format_risk(results['toxicity']['herg_inhibition'])}

Overall: {format_risk(results['toxicity'].get('overall_rating', 'N/A'))}

PHARMACODYNAMICS
IC50: {results['pharmacodynamics']['ic50']:.2f} nM {format_check(results['pharmacodynamics']['ic50'], 1000, False)}
EC50: {results['pharmacodynamics']['ec50']:.2f} nM {format_check(results['pharmacodynamics']['ec50'], 1000, False)}
Binding Affinity: {results['pharmacodynamics']['binding_affinity']:.2f} nM {format_check(results['pharmacodynamics']['binding_affinity'], 500, False)}

Overall: {format_risk(results['pharmacodynamics'].get('overall_rating', 'N/A'))}

ADME PROPERTIES
Absorption:
  Caco-2 Permeability: {results['adme']['absorption']['caco2_permeability']:.2f} nm/s {format_check(results['adme']['absorption']['caco2_permeability'], 20)}
  Oral Absorption: {results['adme']['absorption']['oral_absorption']:.1f}%

Distribution:
  BBB Permeability: {results['adme']['distribution']['bbb_permeability']}
  Plasma Protein Binding: {results['adme']['distribution']['protein_binding']:.1f}%

Metabolism:
  CYP450 Metabolism: {results['adme']['metabolism']['cyp_metabolism']}
  Half-life: {results['adme']['metabolism']['half_life']:.1f} hours

Overall: {format_risk(results['adme'].get('overall_rating', 'N/A'))}

SYNTHETIC ACCESSIBILITY
Synthetic Accessibility Score: {results['stability']['synthetic_accessibility']:.1f} (1-10) {format_check(results['stability']['synthetic_accessibility'], 5, False)}
Estimated Synthetic Steps: {int(results['stability']['synthetic_accessibility'] * 1.2)}

Synthesis Difficulty: {format_risk(results['stability'].get('synthesis_rating', 'N/A'))}

MOLECULAR STABILITY & SHELF LIFE
Water Solubility: {results['stability']['water_solubility']}
Lipid Solubility: {results['stability']['lipid_solubility']}
Shelf Life: {results['stability']['shelf_life']:.1f} months {format_check(results['stability']['shelf_life'], 12)}
Photo-degradation Risk: {format_risk(results['stability']['photodegradation'])}

Overall: {format_risk(results['stability'].get('overall_rating', 'N/A'))}

Overall Drug Candidate Score: {results['overall_assessment']['drug_candidate_score']:.1f}/10

Strengths:
{chr(10).join('- ' + s for s in results['overall_assessment']['strengths'])}

Weaknesses:
{chr(10).join('- ' + w for w in results['overall_assessment']['weaknesses'])}

Final Assessment: {results['overall_assessment']['final_assessment']}"""
    return report

# Update the evaluate_molecule endpoint with better error handling and logging
@app.post("/api/evaluate")
async def evaluate_molecule(request: MoleculeEvalRequest):
    """Evaluate a molecule using ML metrics and generate a detailed report."""
    try:
        logger.info(f"Evaluating molecule with SMILES: {request.smiles}")
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            logger.error(f"Invalid SMILES string: {request.smiles}")
            raise HTTPException(
                status_code=400, 
                detail="Invalid SMILES string - Could not parse molecule"
            )

        # Perform evaluation
        try:
            evaluation_results = predictor.evaluate_molecule(request.smiles)
            logger.info("Molecule evaluation completed successfully")
            
            # Generate report
            report = format_evaluation_report(evaluation_results)
            logger.info("Evaluation report generated")

            # Generate 2D molecule image
            img = Draw.MolToImage(mol)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            
            return JSONResponse(content={
                "metrics": evaluation_results,
                "report": report,
                "smiles": request.smiles,
                "molecule_image": img_str
            })
            
        except Exception as eval_error:
            logger.error(f"Evaluation error: {str(eval_error)}")
            raise HTTPException(
                status_code=400,
                detail=f"Error during molecule evaluation: {str(eval_error)}"
            )
            
    except Exception as e:
        logger.error(f"Error processing evaluation request: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Server error during evaluation: {str(e)}"
        )


if __name__ == "__main__":
    os.makedirs(MODEL_DIR, exist_ok=True)
    if not initialize_models():
        logger.error("Failed to initialize models")
        sys.exit(1)
    uvicorn.run("main:app", host="127.0.0.1", port=5000, reload=True)
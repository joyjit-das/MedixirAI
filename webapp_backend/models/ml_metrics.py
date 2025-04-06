import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, MolSurf, AllChem, rdMolDescriptors
# from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from mordred import Calculator, descriptors
import warnings
warnings.filterwarnings('ignore')
from rdkit.Chem.Crippen import MolLogP
import random  # <-- Add this import

class DrugDiscoveryMetricsPredictor:
    """
    A comprehensive tool to predict various drug discovery metrics from SMILES strings.
    Includes Lipinski's Rule of 5, bioavailability, toxicity, ADME properties, and more.
    """
    
    def __init__(self):
        """Initialize the predictor with necessary models"""
        print("Initializing Drug Discovery Metrics Predictor...")
        # We'll pretend these models are pre-trained and loaded
        # In a real application, you would load actual trained models
        self.bioavailability_model = self._mock_bioavailability_model()
        self.toxicity_models = self._mock_toxicity_models()
        self.adme_models = self._mock_adme_models()
        self.pharmacodynamics_model = self._mock_pharmacodynamics_model()
        # Mordred calculator for additional molecular descriptors
        self.calc = Calculator(descriptors, ignore_3D=True)
        print("Initialization complete!")
        
    def _mock_bioavailability_model(self):
        """Mock bioavailability prediction model"""
        return {
            'oral_bioavailability': RandomForestRegressor(n_estimators=10, random_state=42),
            'first_pass_metabolism': RandomForestRegressor(n_estimators=10, random_state=42),
        }
    
    def _mock_toxicity_models(self):
        """Mock toxicity prediction models"""
        return {
            'hepatotoxicity': RandomForestClassifier(n_estimators=10, random_state=42),
            'nephrotoxicity': RandomForestClassifier(n_estimators=10, random_state=42),
            'carcinogenicity': RandomForestClassifier(n_estimators=10, random_state=42),
            'herg_inhibition': RandomForestClassifier(n_estimators=10, random_state=42),
            'ld50': RandomForestRegressor(n_estimators=10, random_state=42),
        }
    
    def _mock_adme_models(self):
        """Mock ADME property prediction models"""
        return {
            'caco2_permeability': RandomForestRegressor(n_estimators=10, random_state=42),
            'bbb_permeability': RandomForestClassifier(n_estimators=10, random_state=42),
            'plasma_protein_binding': RandomForestRegressor(n_estimators=10, random_state=42),
            'cyp450_inhibition': RandomForestClassifier(n_estimators=10, random_state=42),
            'half_life': RandomForestRegressor(n_estimators=10, random_state=42),
        }
        
    def _mock_pharmacodynamics_model(self):
        """Mock pharmacodynamics prediction models"""
        return {
            'ic50': RandomForestRegressor(n_estimators=10, random_state=42),
            'ec50': RandomForestRegressor(n_estimators=10, random_state=42),
            'kd': RandomForestRegressor(n_estimators=10, random_state=42),
        }
    
    def calculate_lipinski(self, mol):
        """Calculate Lipinski's Rule of Five parameters"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Check how many rules are satisfied
        rules_satisfied = 0
        rules_satisfied += 1 if mw <= 500 else 0
        rules_satisfied += 1 if logp <= 5 else 0
        rules_satisfied += 1 if h_donors <= 5 else 0
        rules_satisfied += 1 if h_acceptors <= 10 else 0
        rules_satisfied += 1 if rotatable_bonds <= 10 else 0
        
        is_drug_like = rules_satisfied >= 4
        
        return {
            'molecular_weight': round(mw, 2),
            'logp': round(logp, 2),
            'h_bond_donors': h_donors,
            'h_bond_acceptors': h_acceptors,
            'rotatable_bonds': rotatable_bonds,
            'rules_satisfied': rules_satisfied,
            'is_drug_like': is_drug_like
        }
    
    def predict_bioavailability(self, mol):
        """Predict bioavailability metrics"""
        # Extract features for prediction
        # In a real model, you would use actual molecular descriptors
        # Here we'll use some basic RDKit descriptors as an example
        features = np.array([
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumRotatableBonds(mol)
        ]).reshape(1, -1)
        
        # Predict using our mock models
        oral_bioavailability = max(0, min(100, 30 + 20 * np.random.randn() + 
                                          0.1 * Descriptors.MolLogP(mol) - 
                                          0.05 * Descriptors.MolWt(mol) / 100))
        
        first_pass = max(0, min(100, 40 + 15 * np.random.randn() + 
                              0.2 * Descriptors.TPSA(mol) / 100))
        
        logd = Descriptors.MolLogP(mol) - 0.5  # Approximation
        
        return {
            'oral_bioavailability': round(oral_bioavailability, 2),
            'first_pass_metabolism': round(first_pass, 2),
            'logd': round(logd, 2),
            'bioavailability_status': 'Good' if oral_bioavailability > 30 else 'Poor'
        }
    
    def predict_toxicity(self, mol):
        """Predict various toxicity parameters"""
        # Extract features for toxicity prediction
        features = np.array([
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.NumAromaticRings(mol),
            Descriptors.NumAliphaticRings(mol)
        ]).reshape(1, -1)
        
        # Simulate toxicity predictions
        # In a real application, these would come from trained models
        tpsa = Descriptors.TPSA(mol)
        logp = Descriptors.MolLogP(mol)
        
        # Higher logP and lower TPSA often correlate with higher hepatotoxicity
        hepatotoxicity_risk = 'High' if logp > 5 or tpsa < 60 else 'Medium' if logp > 3 else 'Low'
        
        # Kidney toxicity often correlates with high molecular weight
        nephrotoxicity_risk = 'High' if Descriptors.MolWt(mol) > 500 and logp > 4 else 'Medium' if Descriptors.MolWt(mol) > 400 else 'Low'
        
        # Carcinogenicity often correlates with specific structural elements
        # This is a simplification
        carcinogenicity_risk = 'Medium' if Descriptors.NumAromaticRings(mol) > 3 else 'Low'
        
        # hERG inhibition correlates with high lipophilicity and presence of basic nitrogen
        n_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N')
        herg_risk = 'High' if logp > 3.7 and n_count > 2 else 'Medium' if logp > 3 else 'Low'
        
        # LD50 approximation (very simplified)
        # Real models would be much more sophisticated
        base_ld50 = 1000 - 100 * logp - Descriptors.MolWt(mol) / 5 + tpsa
        ld50 = max(50, min(2000, base_ld50 + 200 * np.random.randn()))
        
        return {
            'hepatotoxicity': hepatotoxicity_risk,
            'nephrotoxicity': nephrotoxicity_risk,
            'carcinogenicity': carcinogenicity_risk,
            'herg_inhibition': herg_risk,
            'ld50': round(ld50, 2),
            'toxicity_status': 'Safe' if (
                hepatotoxicity_risk == 'Low' and 
                carcinogenicity_risk == 'Low' and 
                ld50 > 200
            ) else 'Caution'
        }
    
    def predict_pharmacodynamics(self, mol):
        """Predict pharmacodynamics parameters (IC50, EC50, Kd)"""
        # These would normally be predicted by specialized ML models
        # Here we're simulating based on molecular properties
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # More complex molecules with balanced properties often have better binding
        complexity = Descriptors.BalabanJ(mol) if Descriptors.BalabanJ(mol) > 0 else 1
        
        # IC50 simulation (nanomolar range)
        base_ic50 = 500 - 0.5 * tpsa + 100 * logp / complexity
        ic50 = max(0.1, base_ic50 + 100 * np.random.randn())
        
        # EC50 simulation (nanomolar range)
        ec50 = ic50 * (0.5 + np.random.rand())
        
        # Kd simulation (nanomolar range)
        kd = ic50 * (0.3 + 0.4 * np.random.rand())
        
        return {
            'ic50': round(ic50, 2),
            'ec50': round(ec50, 2),
            'kd': round(kd, 2),
            'binding_quality': 'Strong' if ic50 < 1000 and kd < 10 else 'Moderate' if ic50 < 10000 else 'Weak'
        }
    
    def predict_adme(self, mol):
        """Predict ADME properties"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        
        # Caco-2 permeability prediction (nm/s)
        # Higher logP and lower TPSA generally correlate with higher permeability
        caco2 = max(0, 20 - 0.1 * tpsa + 5 * logp - 0.5 * h_donors - 0.5 * h_acceptors + 5 * np.random.randn())
        
        # Blood-Brain Barrier permeability
        # Simple rule: if MW < 400, logP between 1-4, TPSA < 90, likely to cross BBB
        bbb_permeability = 'High' if (mw < 400 and 1 < logp < 4 and tpsa < 90) else 'Low'
        
        # Plasma Protein Binding
        # Higher logP often correlates with higher binding
        ppb = min(99, max(10, 50 + 10 * logp + 10 * np.random.randn()))
        
        # CYP450 metabolism prediction
        # This is highly complex, but higher logP molecules often have higher metabolism
        cyp_metabolism = 'High' if logp > 4 else 'Medium' if logp > 2 else 'Low'
        
        # Half-life prediction (hours)
        half_life = max(0.5, min(24, 6 + 3 * np.random.randn() - logp + tpsa / 100))
        
        return {
            'caco2_permeability': round(caco2, 2),
            'bbb_permeability': bbb_permeability,
            'plasma_protein_binding': round(ppb, 2),
            'cyp450_metabolism': cyp_metabolism,
            'half_life': round(half_life, 2),
            'adme_quality': 'Good' if (caco2 > 10 and 1 <= logp <= 3 and cyp_metabolism != 'High' and 3 <= half_life <= 12) else 'Suboptimal'
        }
    
    def predict_synthetic_accessibility(self, mol):
        """Predict synthetic accessibility and related metrics"""
        # RDKit's synthetic accessibility score
        try:
            sas = round(Chem.rdMolDescriptors.CalcSyntheticAccessibilityScore(mol), 2)
        except:
            sas = 5.0  # Default if calculation fails
        
        # Estimate number of synthetic steps (very approximate)
        # Based on molecule complexity and ring count
        complexity = Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol) + Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)
        ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
        estimated_steps = max(1, round(1 + sas / 2 + ring_count / 2 + complexity))
        
        return {
            'synthetic_accessibility_score': sas,
            'estimated_synthetic_steps': estimated_steps,
            'synthesis_difficulty': 'Easy' if sas < 5 else 'Moderate' if sas < 7 else 'Difficult'
        }
    
    def predict_stability(self, mol):
        """Predict molecular stability and shelf life"""
        # This is a complex property to predict, so we're using some approximations
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # More rotatable bonds and extreme logP values can lead to lower stability
        stability_factor = 5 - abs(logp - 2) / 2 - rotatable_bonds / 10
        
        # Water solubility prediction
        # Lower logP correlates with higher water solubility
        water_solubility = 'High' if logp < 0 else 'Moderate' if logp < 3 else 'Low'
        
        # Lipid solubility is roughly the opposite
        lipid_solubility = 'Low' if logp < 0 else 'Moderate' if logp < 3 else 'High'
        
        # Degradation half-life (months)
        # This is affected by many factors, this is a very rough estimate
        if stability_factor > 3:
            shelf_life = 18 + 6 * np.random.rand()  # 18-24 months
        elif stability_factor > 1:
            shelf_life = 12 + 6 * np.random.rand()  # 12-18 months
        else:
            shelf_life = 6 + 6 * np.random.rand()   # 6-12 months
        
        # Photo-degradation risk
        # Higher for molecules with certain functional groups (simplified)
        photo_degradation = 'Low'
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ['N', 'O', 'S'] and atom.GetDegree() == 1:
                photo_degradation = 'Moderate'
                break
        
        return {
            'water_solubility': water_solubility,
            'lipid_solubility': lipid_solubility,
            'shelf_life_months': round(shelf_life, 1),
            'photo_degradation_risk': photo_degradation,
            'stability_status': 'Stable' if shelf_life > 12 else 'Moderate' if shelf_life > 6 else 'Unstable'
        }
    
    def predict_all_metrics(self, smiles):
        """Predict all drug discovery metrics for a given SMILES string"""
        try:
            # Convert SMILES to RDKit molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": "Invalid SMILES string"}
            
            # Calculate all metrics
            results = {
                "smiles": smiles,
                "molecule_name": "Compound",  # Default name
                "structure_info": {
                    "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
                    "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
                    "heavy_atom_count": mol.GetNumHeavyAtoms(),
                    "ring_count": Chem.rdMolDescriptors.CalcNumRings(mol)
                },
            }
            
            # Calculate each category of metrics
            results["drug_likeness"] = self.calculate_lipinski(mol)
            results["bioavailability"] = self.predict_bioavailability(mol)
            results["toxicity"] = self.predict_toxicity(mol)
            results["pharmacodynamics"] = self.predict_pharmacodynamics(mol)
            results["adme_properties"] = self.predict_adme(mol)
            results["synthetic_accessibility"] = self.predict_synthetic_accessibility(mol)
            results["stability"] = self.predict_stability(mol)
            
            # Calculate QED (Quantitative Estimate of Drug-likeness)
            results["qed"] = round(QED.qed(mol), 3)
            
            return results
            
        except Exception as e:
            return {"error": f"Error processing SMILES: {str(e)}"}
    
    def generate_report(self, results):
        """Generate a formatted text report from the prediction results without emojis."""
        if "error" in results:
            return f"Error: {results['error']}"
        
        report = f"""
╔══════════════════════════════════════════════════════════════╗
║                 DRUG DISCOVERY METRICS REPORT                 ║
╚══════════════════════════════════════════════════════════════╝

COMPOUND INFORMATION
-------------------
SMILES: {results['smiles']}
Formula: {results['structure_info']['formula']}
Exact Mass: {results['structure_info']['exact_mass']} Da
Heavy Atoms: {results['structure_info']['heavy_atom_count']}
Ring Count: {results['structure_info']['ring_count']}
Overall QED Score: {results['qed']}

DRUG-LIKENESS (Lipinski's Rule of Five)
-------------------------------------------
Molecular Weight: {results['drug_likeness']['molecular_weight']} Da {"(Pass)" if results['drug_likeness']['molecular_weight'] <= 500 else "(Fail)"}
LogP: {results['drug_likeness']['logp']} {"(Pass)" if results['drug_likeness']['logp'] <= 5 else "(Fail)"}
H-Bond Donors: {results['drug_likeness']['h_bond_donors']} {"(Pass)" if results['drug_likeness']['h_bond_donors'] <= 5 else "(Fail)"}
H-Bond Acceptors: {results['drug_likeness']['h_bond_acceptors']} {"(Pass)" if results['drug_likeness']['h_bond_acceptors'] <= 10 else "(Fail)"}
Rotatable Bonds: {results['drug_likeness']['rotatable_bonds']} {"(Pass)" if results['drug_likeness']['rotatable_bonds'] <= 10 else "(Fail)"}

Rules Satisfied: {results['drug_likeness']['rules_satisfied']}/5
Overall: {"DRUG-LIKE" if results['drug_likeness']['is_drug_like'] else "NOT DRUG-LIKE"}

BIOAVAILABILITY
------------------
Oral Bioavailability: {results['bioavailability']['oral_bioavailability']}% {"(Pass)" if results['bioavailability']['oral_bioavailability'] > 30 else "(Fail)"}
First-Pass Metabolism: {results['bioavailability']['first_pass_metabolism']}%
LogD: {results['bioavailability']['logd']}

Overall: {"GOOD" if results['bioavailability']['bioavailability_status'] == 'Good' else "POOR"}

TOXICITY PREDICTION
----------------------
Hepatotoxicity Risk: {results['toxicity']['hepatotoxicity']} {"(Low)" if results['toxicity']['hepatotoxicity'] == "Low" else "(Medium)" if results['toxicity']['hepatotoxicity'] == "Medium" else "(High)"}
Nephrotoxicity Risk: {results['toxicity']['nephrotoxicity']} {"(Low)" if results['toxicity']['nephrotoxicity'] == "Low" else "(Medium)" if results['toxicity']['nephrotoxicity'] == "Medium" else "(High)"}
Carcinogenicity Risk: {results['toxicity']['carcinogenicity']} {"(Low)" if results['toxicity']['carcinogenicity'] == "Low" else "(Medium)" if results['toxicity']['carcinogenicity'] == "Medium" else "(High)"}
hERG Inhibition Risk: {results['toxicity']['herg_inhibition']} {"(Low)" if results['toxicity']['herg_inhibition'] == "Low" else "(Medium)" if results['toxicity']['herg_inhibition'] == "Medium" else "(High)"}
LD50: {results['toxicity']['ld50']} mg/kg {"(Pass)" if results['toxicity']['ld50'] > 200 else "(Fail)"}

Overall: {"SAFE" if results['toxicity']['toxicity_status'] == 'Safe' else "CAUTION"}

PHARMACODYNAMICS
-------------------
IC50: {results['pharmacodynamics']['ic50']} nM {"(Pass)" if results['pharmacodynamics']['ic50'] < 1000 else "(Fail)"}
EC50: {results['pharmacodynamics']['ec50']} nM {"(Pass)" if results['pharmacodynamics']['ec50'] < 1000 else "(Fail)"}
Kd: {results['pharmacodynamics']['kd']} nM {"(Pass)" if results['pharmacodynamics']['kd'] < 10 else "(Fail)"}

Overall Binding: {"STRONG" if results['pharmacodynamics']['binding_quality'] == 'Strong' else "MODERATE" if results['pharmacodynamics']['binding_quality'] == 'Moderate' else "WEAK"}

ADME PROPERTIES
------------------
Absorption:
  Caco-2 Permeability: {results['adme_properties']['caco2_permeability']} nm/s {"(Pass)" if results['adme_properties']['caco2_permeability'] > 10 else "(Fail)"}
  LogP: {results['drug_likeness']['logp']} {"(Pass)" if 1 <= results['drug_likeness']['logp'] <= 3 else "(Fail)"}

Distribution:
  BBB Permeability: {results['adme_properties']['bbb_permeability']} 
  Plasma Protein Binding: {results['adme_properties']['plasma_protein_binding']}%

Metabolism:
  CYP450 Metabolism: {results['adme_properties']['cyp450_metabolism']} {"(Pass)" if results['adme_properties']['cyp450_metabolism'] != 'High' else "(Fail)"}

Excretion:
  Half-life: {results['adme_properties']['half_life']} hours {"(Pass)" if 3 <= results['adme_properties']['half_life'] <= 12 else "(Fail)"}

Overall: {"GOOD" if results['adme_properties']['adme_quality'] == 'Good' else "SUBOPTIMAL"}

SYNTHETIC ACCESSIBILITY
--------------------------
Synthetic Accessibility Score: {results['synthetic_accessibility']['synthetic_accessibility_score']} (1-10) {"(Pass)" if results['synthetic_accessibility']['synthetic_accessibility_score'] < 5 else "(Fail)"}
Estimated Synthetic Steps: {results['synthetic_accessibility']['estimated_synthetic_steps']}

Synthesis Difficulty: {"EASY" if results['synthetic_accessibility']['synthesis_difficulty'] == 'Easy' else "MODERATE" if results['synthetic_accessibility']['synthesis_difficulty'] == 'Moderate' else "DIFFICULT"}

MOLECULAR STABILITY & SHELF LIFE
----------------------------------
Water Solubility: {results['stability']['water_solubility']}
Lipid Solubility: {results['stability']['lipid_solubility']}
Shelf Life: {results['stability']['shelf_life_months']} months {"(Pass)" if results['stability']['shelf_life_months'] > 12 else "(Fail)"}
Photo-degradation Risk: {results['stability']['photo_degradation_risk']} {"(Pass)" if results['stability']['photo_degradation_risk'] == 'Low' else "(Fail)"}

Overall: {"STABLE" if results['stability']['stability_status'] == 'Stable' else "MODERATE" if results['stability']['stability_status'] == 'Moderate' else "UNSTABLE"}


"""
        
        # Create a summary scoring system
        drug_like_score = 1 if results['drug_likeness']['is_drug_like'] else 0
        bioavail_score = 1 if results['bioavailability']['bioavailability_status'] == 'Good' else 0
        toxicity_score = 1 if results['toxicity']['toxicity_status'] == 'Safe' else 0.5 if results['toxicity']['toxicity_status'] == 'Caution' else 0
        pd_score = 1 if results['pharmacodynamics']['binding_quality'] == 'Strong' else 0.5 if results['pharmacodynamics']['binding_quality'] == 'Moderate' else 0
        adme_score = 1 if results['adme_properties']['adme_quality'] == 'Good' else 0.5
        synth_score = 1 if results['synthetic_accessibility']['synthesis_difficulty'] == 'Easy' else 0.5 if results['synthetic_accessibility']['synthesis_difficulty'] == 'Moderate' else 0
        stability_score = 1 if results['stability']['stability_status'] == 'Stable' else 0.5 if results['stability']['stability_status'] == 'Moderate' else 0
        
        total_score = (drug_like_score + bioavail_score + toxicity_score + pd_score + adme_score + synth_score + stability_score) / 7 * 10
        
        # Add summary assessment
        report += f"""
Overall Drug Candidate Score: {total_score:.1f}/10

Strengths:
"""
        if drug_like_score == 1:
            report += "- Follows Lipinski's Rule of Five criteria\n"
        if bioavail_score == 1:
            report += "- Good bioavailability profile\n"
        if toxicity_score == 1:
            report += "- Low toxicity risks\n"
        if pd_score == 1:
            report += "- Strong binding affinity to targets\n"
        if adme_score == 1:
            report += "- Favorable ADME properties\n"
        if synth_score == 1:
            report += "- Easy to synthesize\n"
        if stability_score == 1:
            report += "- Good stability profile\n"
        
        report += "\nWeaknesses:\n"
        if drug_like_score == 0:
            report += "- Does not meet Lipinski's Rule of Five criteria\n"
        if bioavail_score == 0:
            report += "- Poor bioavailability\n"
        if toxicity_score < 1:
            report += "- Potential toxicity concerns\n"
        if pd_score < 1:
            report += "- Suboptimal target binding\n"
        if adme_score < 1:
            report += "- Suboptimal ADME properties\n"
        if synth_score < 1:
            report += "- Challenging synthesis\n"
        if stability_score < 1:
            report += "- Stability concerns\n"
        
        report += f"""
Final Assessment: This compound is {"PROMISING" if total_score >= 7 else "MODERATE" if total_score >= 5 else "NOT PROMISING"} as a drug candidate.
"""
        
        return report

    def evaluate_molecule(self, smiles):
        """Main evaluation method"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string")

            # Calculate basic properties
            mw = Descriptors.ExactMolWt(mol)
            logp = MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            rotatable = Descriptors.NumRotatableBonds(mol)
            qed = QED.default(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            rings = rdMolDescriptors.CalcNumRings(mol)

            # Structure results
            results = {
                "compound_info": {
                    "smiles": smiles,
                    "formula": formula,
                    "exact_mass": mw,
                    "heavy_atoms": mol.GetNumHeavyAtoms(),
                    "ring_count": rings,
                    "qed_score": float(qed)
                },
                "drug_likeness": {
                    "molecular_weight": mw,
                    "logp": logp,
                    "hbd": hbd,
                    "hba": hba,
                    "rotatable_bonds": rotatable,
                    "lipinski_violations": sum([
                        mw > 500,
                        logp > 5,
                        hbd > 5,
                        hba > 10
                    ]),
                    "is_drug_like": True  # Will be updated
                },
                "bioavailability": {
                    "oral_bioavailability": max(0, min(100, 65 - abs(logp - 2) * 10 - tpsa/20)),
                    "first_pass_metabolism": min(95, max(5, 50 + logp * 5)),
                    "logd": logp - 0.5,
                    "overall_rating": "Good"  # Will be updated
                },
                "toxicity": {
                    "hepatotoxicity": "Low" if logp <= 3 else "Medium" if logp <= 5 else "High",
                    "nephrotoxicity": "Low" if mw <= 400 else "Medium" if mw <= 500 else "High",
                    "carcinogenicity": "Low" if rings <= 3 else "Medium" if rings <= 4 else "High",
                    "herg_inhibition": "Low" if logp <= 3 else "Medium" if logp <= 4 else "High",
                    "overall_rating": "Low"  # Will be updated
                },
                "pharmacodynamics": {
                    "ic50": 500 + random.uniform(-200, 200),
                    "ec50": 450 + random.uniform(-150, 150),
                    "binding_affinity": 300 + random.uniform(-100, 100),
                    "overall_rating": "Moderate"
                },
                "adme": {
                    "absorption": {
                        "caco2_permeability": max(0, 20 - tpsa/10 + logp*2),
                        "oral_absorption": min(100, max(0, 80 - tpsa/5))
                    },
                    "distribution": {
                        "bbb_permeability": "Low" if tpsa > 90 else "High",
                        "protein_binding": min(99.9, max(10, 50 + logp * 10))
                    },
                    "metabolism": {
                        "cyp_metabolism": "High" if logp > 4 else "Medium" if logp > 2 else "Low",
                        "half_life": max(0.5, min(24, 12 - abs(logp-2)*2))
                    }
                },
                "stability": {
                    "water_solubility": "Low" if logp > 3 else "Medium" if logp > 0 else "High",
                    "lipid_solubility": "High" if logp > 2 else "Medium" if logp > 0 else "Low",
                    "shelf_life": max(6, min(24, 18 - abs(logp-2)*2)),
                    "photodegradation": "Low",
                    "synthetic_accessibility": min(10, max(1, 5 + logp/2)),
                    "synthesis_rating": "Moderate"
                }
            }

            # Update derived values
            results["drug_likeness"]["is_drug_like"] = results["drug_likeness"]["lipinski_violations"] <= 1
            results["bioavailability"]["overall_rating"] = "Good" if results["bioavailability"]["oral_bioavailability"] > 30 else "Poor"
            results["toxicity"]["overall_rating"] = "Low" if all(x == "Low" for x in [results["toxicity"]["hepatotoxicity"], results["toxicity"]["nephrotoxicity"]]) else "High"
            
            # Calculate overall assessment
            strengths = []
            if results["drug_likeness"]["is_drug_like"]:
                strengths.append("Follows Lipinski's Rule of Five criteria")
            if results["bioavailability"]["overall_rating"] == "Good":
                strengths.append("Good bioavailability profile")
            if results["toxicity"]["overall_rating"] == "Low":
                strengths.append("Low toxicity risk")

            weaknesses = []
            if not results["drug_likeness"]["is_drug_like"]:
                weaknesses.append("Does not meet Lipinski's Rule of Five criteria")
            if results["bioavailability"]["overall_rating"] == "Poor":
                weaknesses.append("Poor bioavailability")
            if results["toxicity"]["overall_rating"] == "High":
                weaknesses.append("High toxicity risk")

            score = 7.0 + len(strengths) - len(weaknesses)
            score = max(1, min(10, score))

            results["overall_assessment"] = {
                "drug_candidate_score": score,
                "strengths": strengths,
                "weaknesses": weaknesses,
                "final_assessment": "PROMISING" if score >= 7 else "MODERATE" if score >= 5 else "NOT PROMISING"
            }

            return results

        except Exception as e:
            raise ValueError(f"Error evaluating molecule: {str(e)}")

def main():
    """Main function to demonstrate the drug metrics predictor"""
    predictor = DrugDiscoveryMetricsPredictor()
    
    print("\nDrug Discovery Metrics Predictor")
    print("===============================")
    print("This tool predicts various drug discovery metrics from SMILES strings.")
    
    while True:
        print("\nEnter a SMILES string (or 'q' to quit):")
        smiles = input("> ")
        
        if smiles.lower() == 'q':
            break
        
        # Example SMILES if none provided
        if not smiles:
            smiles = "CCO"  # Ethanol
            print(f"Using example SMILES: {smiles}")
        
        # Predict metrics
        results = predictor.predict_all_metrics(smiles)
        
        # Display report
        report = predictor.generate_report(results)
        print(report)

    
    print("Thank you for using the Drug Discovery Metrics Predictor!")

if __name__ == "__main__":
    main()

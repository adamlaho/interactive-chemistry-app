import streamlit as st
import random
import pandas as pd
import time

# RDKit and py3Dmol imports for molecular processing and visualization
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

# 🏆 Streamlit Title
st.title("🧪 Interactive Chemistry App")

# ===================== 1. Generate a Random Molecule (with Lewis Structure) =====================
st.header("🔬 Generate a Random Molecule (with Lewis Structure)")
st.write("💡 Select a set of atoms to create a random molecule:")
atom_choices = ["H", "O", "N", "C", "F", "Cl", "Br", "S"]
selected_atoms = st.multiselect("Choose atoms for the molecule:", atom_choices, default=["C", "H", "O"])

def generate_random_molecule(selected_atoms, num_atoms=None):
    """
    Generate a molecule from a random chain of atoms.
    """
    smiles_list = {
        "H": "[H]",
        "O": "O",
        "N": "N",
        "C": "C",
        "F": "F",
        "Cl": "Cl",
        "Br": "Br",
        "S": "S"
    }
    # If num_atoms not provided, default to a small molecule (2-6 atoms)
    if num_atoms is None:
        num_atoms = random.randint(2, 6)
    selected_smiles = random.choices([smiles_list[a] for a in selected_atoms], k=num_atoms)
    smiles = "".join(selected_smiles)
    mol = Chem.MolFromSmiles(smiles)
    return mol, smiles

if st.button("Generate Molecule (2D)"):
    mol, smiles = generate_random_molecule(selected_atoms)
    if mol:
        st.image(Draw.MolToImage(mol), caption=f"Lewis Structure of: {smiles}")
    else:
        st.error("Could not generate a valid molecule. Try different atoms!")

# ===================== 2. Interactive Periodic Table =====================
st.header("📊 Interactive Periodic Table")
periodic_table = pd.DataFrame({
    "Atomic Number": list(range(1, 19)),
    "Element": ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"],
    "Group": ["Non-metal", "Noble Gas", "Alkali Metal", "Alkaline Earth Metal", "Metalloid", "Non-metal",
              "Non-metal", "Non-metal", "Halogen", "Noble Gas", "Alkali Metal", "Alkaline Earth Metal",
              "Metal", "Metalloid", "Non-metal", "Non-metal", "Halogen", "Noble Gas"]
})
num_elements = st.slider("Select number of elements to display:", 1, 18, 10)
filtered_table = periodic_table.head(num_elements)
st.dataframe(filtered_table)

# ===================== 3. Interactive 3D Molecule Viewer (Big Molecules in Ball-and-Stick) =====================
st.header("🔍 Interactive 3D Molecule Viewer (Big Molecules)")
st.write("Select the number of atoms for your 3D molecule:")
num_atoms_3d = st.slider("Number of atoms:", min_value=1, max_value=15, value=12)

if st.button("Generate 3D Molecule"):
    st.info("Searching for a valid 3D molecule.")
    # Keep trying until a valid molecule is found
    attempt_count = 0
    while True:
        attempt_count += 1
        mol, smiles = generate_random_molecule(selected_atoms, num_atoms=num_atoms_3d)
        if mol is None:
            continue  # Try again if no molecule was generated
        mol = Chem.AddHs(mol)
        status = AllChem.EmbedMolecule(mol)
        if status == 0:
            # Found a valid molecule
            break
        # Optional: Sleep a short time to avoid rapid-fire looping
        time.sleep(0.1)
    
    #st.success(f"Valid molecule found after {attempt_count} attempts!")
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    
    # Create a 3D viewer using py3Dmol with ball-and-stick style.
    view = py3Dmol.view(width=600, height=600)
    view.addModel(pdb_block, 'pdb')
    view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    view_html = view._make_html()
    st.components.v1.html(view_html, height=450)

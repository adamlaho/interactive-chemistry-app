import streamlit as st
import random
import pandas as pd
import time
import pubchempy as pcp

# RDKit and py3Dmol imports for molecular processing and visualization
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, rdMolDescriptors
import py3Dmol
import streamlit.components.v1 as components

# Streamlit title
st.title("🧪 Interactive Chemistry App")

# ===================== 1. Generate a Random Molecule (with Lewis Structure) =====================
st.header("🔬 Generate a Random Molecule (with Lewis Structure)")
st.write("💡 Select a set of atoms to create a random molecule (up to 50 attempts).")

atom_choices = ["H", "O", "N", "C", "F", "Cl", "Br", "S"]
selected_atoms = st.multiselect("Choose atoms for the molecule:", atom_choices, default=["C", "H", "O"])

def generate_random_smiles(selected_atoms, num_atoms=5):
    """
    Generate a random SMILES string by concatenating a given number of atoms
    randomly selected from the provided list.
    """
    if not selected_atoms:
        return ""
    smiles = ""
    for i in range(num_atoms):
        atom = random.choice(selected_atoms)
        smiles += atom
    return smiles

if selected_atoms:
    # Attempt multiple times to get a valid molecule
    max_attempts = 50
    valid_mol = None
    random_smiles_final = None

    for attempt in range(max_attempts):
        candidate_smiles = generate_random_smiles(selected_atoms, num_atoms=5)
        candidate_mol = Chem.MolFromSmiles(candidate_smiles)
        if candidate_mol:
            valid_mol = candidate_mol
            random_smiles_final = candidate_smiles
            break

    if valid_mol:
        st.write(f"**Generated SMILES** (valid after {attempt+1} attempt(s)): {random_smiles_final}")
        st.subheader("Lewis Structure")
        img = Draw.MolToImage(valid_mol, size=(300, 300))
        st.image(img)
    else:
        st.error(f"Could not generate a valid molecule after {max_attempts} attempts.")

# ===================== 2. Biological Properties Evaluation =====================
st.header("🧬 Biological Properties Evaluation")
st.write("Select an amino acid to visualize its 3D structure and view computed properties:")

# Dictionary of amino acids with canonical SMILES representations
aa_smiles = {
    "Alanine": "N[C@@H](C)C(=O)O",
    "Glycine": "NCC(=O)O",
    "Cysteine": "N[C@@H](CS)C(=O)O",
    "Serine": "N[C@@H](CO)C(=O)O",
    "Valine": "N[C@@H](C(C)C)C(=O)O"
}

aa_choice = st.selectbox("Select an amino acid:", list(aa_smiles.keys()))

if aa_choice:
    smiles = aa_smiles[aa_choice]
    st.write(f"**SMILES:** `{smiles}`")

    def generate_3d_structure(smiles_str):
        """
        Generate a 3D molecular structure:
         - Create molecule from SMILES
         - Add explicit hydrogens
         - Embed in 3D and optimize geometry using UFF
        """
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            st.error("Error generating molecule from SMILES.")
            return None
        mol = Chem.AddHs(mol)
        status = AllChem.EmbedMolecule(mol)
        if status != 0:
            st.error("3D embedding failed.")
            return None
        AllChem.UFFOptimizeMolecule(mol)
        return mol

    mol_aa = generate_3d_structure(smiles)
    if mol_aa:
        def show_molecule(mol):
            """
            Create a 3D visualization using py3Dmol in a unified ball-and-stick style,
            with the correct 'sdf' model type for a MolToMolBlock output.
            """
            mb = Chem.MolToMolBlock(mol)
            view = py3Dmol.view(width=400, height=300)
            view.addModel(mb, 'sdf')
            # Single style for all atoms: ball-and-stick with element-based color scheme
            view.setStyle({},
                {
                    "stick": {
                        "colorscheme": "element",
                        "radius": 0.15
                    },
                    "sphere": {
                        "colorscheme": "element",
                        "scale": 0.3
                    }
                }
            )
            view.zoomTo()
            return view

        view = show_molecule(mol_aa)
        st.subheader("3D Structure")
        components.html(view._make_html(), height=320)

        # Compute and display molecular properties
        properties = {
            "Molecular Weight": round(Descriptors.MolWt(mol_aa), 2),
            "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol_aa),
            "Topological Polar Surface Area": round(Descriptors.TPSA(mol_aa), 2),
            "Hydrogen Bond Donors": Descriptors.NumHDonors(mol_aa),
            "Hydrogen Bond Acceptors": Descriptors.NumHAcceptors(mol_aa),
        }
        st.subheader("Computed Properties")
        for prop, value in properties.items():
            st.write(f"**{prop}:** {value}")

# ===================== 3. Interactive Periodic Table =====================
st.header("📊 Interactive Periodic Table")
periodic_table = pd.DataFrame({
    "Atomic Number": list(range(1, 19)),
    "Element": ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"],
    "Group": ["Non-metal", "Noble Gas", "Alkali Metal", "Alkaline Earth Metal",
              "Metalloid", "Non-metal", "Non-metal", "Non-metal", "Halogen",
              "Noble Gas", "Alkali Metal", "Alkaline Earth Metal", "Metal",
              "Metalloid", "Non-metal", "Non-metal", "Halogen", "Noble Gas"]
})
num_elements = st.slider("Select number of elements to display:", 1, 18, 10)
filtered_table = periodic_table.head(num_elements)
st.dataframe(filtered_table)

# ===================== 4. Interactive 3D Molecule Viewer (From PubChem Database) =====================
st.header("🔍 Interactive 3D Molecule Viewer (PubChem)")
st.write("Search PubChem by a compound name, then display its 3D structure in ball-and-stick style.")

compound_name = st.text_input("Enter a compound name to fetch from PubChem:", "aspirin")

if st.button("Fetch PubChem Molecule"):
    try:
        st.info(f"Querying PubChem for '{compound_name}'...")
        # Search PubChem by name
        compounds = pcp.get_compounds(compound_name, "name")

        if not compounds:
            st.warning("No compounds found. Please try a different query.")
        else:
            chosen_compound = compounds[0]  # Take the first match
            st.success(f"Found molecule: CID {chosen_compound.cid}")
            st.write(f"SMILES: {chosen_compound.canonical_smiles}")

            # Convert the SMILES to RDKit, add explicit hydrogens
            mol_pub = Chem.MolFromSmiles(chosen_compound.canonical_smiles)
            if mol_pub is None:
                st.error("Could not parse the SMILES from PubChem.")
            else:
                mol_pub = Chem.AddHs(mol_pub)
                status = AllChem.EmbedMolecule(mol_pub)
                if status != 0:
                    st.error("3D embedding failed for this molecule.")
                else:
                    AllChem.UFFOptimizeMolecule(mol_pub)
                    # Create an SDF block for Py3Dmol
                    sdf_block = Chem.MolToMolBlock(mol_pub)

                    # Create a Py3Dmol viewer
                    view = py3Dmol.view(width=600, height=600)
                    view.addModel(sdf_block, 'sdf')
                    # Single style for all atoms: ball-and-stick with an element-based color scheme
                    view.setStyle({},
                        {
                            "stick": {
                                "colorscheme": "element",
                                "radius": 0.15
                            },
                            "sphere": {
                                "colorscheme": "element",
                                "scale": 0.3
                            }
                        }
                    )
                    view.zoomTo()

                    st.subheader("3D Structure (PubChem)")
                    components.html(view._make_html(), height=450)
    
                    # Print out some properties of the PubChem molecule
                    pub_props = {
                        "Molecular Weight": round(Descriptors.MolWt(mol_pub), 2),
                        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol_pub),
                        "Topological Polar Surface Area": round(Descriptors.TPSA(mol_pub), 2),
                        "Hydrogen Bond Donors": Descriptors.NumHDonors(mol_pub),
                        "Hydrogen Bond Acceptors": Descriptors.NumHAcceptors(mol_pub),
                    }
                    st.subheader("Computed Properties (PubChem Molecule)")
                    for prop, value in pub_props.items():
                        st.write(f"**{prop}:** {value}")
        
    except pcp.PubChemHTTPError as e:
        st.error(f"PubChem query failed: {e}")
    except Exception as ex:
        st.error(f"An unexpected error occurred: {ex}")

"""Utilities to compute energy-storage-related properties with xTB"""
from io import StringIO

from xtb.ase.calculator import XTB
from ase.optimize import QuasiNewton
from ase.io import read
from ase import units
from rdkit.Chem import AllChem
from rdkit import Chem

from cfree.descriptors import saturate_molecule, count_h2_difference

# Energy of H2 in isolation
_h2_energy = -26.74024588370831


# Adapted from ExaMol: https://github.com/exalearn/ExaMol/blob/main/examol/simulate/initialize.py
def generate_xyz(smiles: str) -> str:
    """Generate the XYZ coordinates for a molecule

    Args:
        smiles: SMILES of molecule
    Returns:
        XYZ coordinates for the molecule
    """

    # Generate 3D coordinates for the molecule
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)

    # Save geometry as 3D coordinates
    xyz = f"{mol.GetNumAtoms()}\n"
    xyz += smiles + "\n"
    conf = mol.GetConformer()
    for i, a in enumerate(mol.GetAtoms()):
        s = a.GetSymbol()
        c = conf.GetAtomPosition(i)
        xyz += f"{s} {c[0]} {c[1]} {c[2]}\n"

    return xyz


def compute_storage_energy(smiles: str) -> float:
    """Compute the storage energy for a molecule

    Compute the relaxed energy of the current state and fully-hydrated state

    Args:
        smiles: SMILES of the molecule in its dehydrated state
    Returns:
        Storage energy in kJ/mol
    """

    # Get the XYZ of the fully-hydrated and current state
    hyd_smiles = saturate_molecule(smiles)
    dehyd_xyz = generate_xyz(smiles)
    hyd_xyz = generate_xyz(hyd_smiles)

    # Make a function to compute the energy of the relaxed molecule
    calc = XTB()

    def _get_energy(xyz: str) -> float:
        """Get the total energy of the relaxed structure"""

        atoms = read(StringIO(xyz), format='xyz')
        atoms.calc = calc
        opt = QuasiNewton(atoms, logfile=None)
        opt.run(fmax=0.01)
        return atoms.get_total_energy()

    # Compute the total energy of the two states
    dehyd_energy = _get_energy(dehyd_xyz)
    hyd_energy = _get_energy(hyd_xyz)

    # Return the energy difference
    h2_count = count_h2_difference(smiles, hyd_smiles)
    rxn_energy = (dehyd_energy + h2_count * _h2_energy) - hyd_energy
    rxn_energy /= h2_count
    rxn_energy /= units.kJ / units.mol
    return rxn_energy

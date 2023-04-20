"""Tools for computing descriptors related to hydrogen storage capability"""
from typing import Optional

from rdkit.Chem import Descriptors
from rdkit import Chem


def saturate_molecule(molecule_smiles: str) -> str:
    """Generate a fully-saturated version of a molecule

    Args:
        molecule_smiles: SMILES string of molecule to be saturated
    Returns:
        SMILES string of saturated molecule
    """
    # Code originally by Hassan Harb
    new_smiles = ''
    for i in molecule_smiles:
        if i == '=' or i == '#':  # Skip double and triple bonds
            continue
        elif i in {'c', 'n', 'o', 's'}:
            i = i.upper()
            new_smiles = new_smiles + i
        else:
            new_smiles = new_smiles + i
    return new_smiles


def compute_wth2(smiles_dehydrogenated: str, smiles_hydrogenated: Optional[str] = None) -> float:
    """Compute the wt% of hydrogen that is stored by a molecule

    Args:
        smiles_dehydrogenated: Dehydrogenated form of the molecule
        smiles_hydrogenated: Hydrogenated form. If `None`, we will compute the fully-hydrogenated form
    Returns:
        wt% of hydrogen stored by molecule
    """
    if smiles_hydrogenated is None:
        smiles_hydrogenated = saturate_molecule(smiles_dehydrogenated)
    mw_hydrogenated = Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles_hydrogenated))
    mw_dehydrogenated = Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles_dehydrogenated))
    return ((mw_hydrogenated - mw_dehydrogenated) / mw_hydrogenated) * 100


def count_h2_difference(smiles_dehydrogenated: str, smiles_hydrogenated: Optional[str] = None) -> int:
    """Count the H2 difference between hydrogenated and unhydrogenated state

    Args:
        smiles_dehydrogenated: Dehydrogenated form of the molecule
        smiles_hydrogenated: Hydrogenated form. If `None`, we will compute the fully-hydrogenated form
    Returns:
        Number of H2s between each
    """

    if smiles_hydrogenated is None:
        smiles_hydrogenated = saturate_molecule(smiles_dehydrogenated)

    # Count H2s
    def _count(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        return len([x for x in mol.GetAtoms() if x.GetAtomicNum() == 1])

    return (_count(smiles_hydrogenated) - _count(smiles_dehydrogenated)) // 2

from cfree.store import MoleculeRecord

from pytest import raises


def test_initialization():
    doc = MoleculeRecord.from_identifier(smiles='C')
    doc_2 = MoleculeRecord.from_identifier(inchi=doc.identifier.inchi)
    assert doc.key == doc_2.key
    assert doc.key is not None

    with raises(ValueError):
        MoleculeRecord.from_identifier(inchi='Nah')

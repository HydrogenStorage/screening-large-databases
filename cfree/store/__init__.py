from typing import Optional

from mongoengine import Document, DynamicEmbeddedDocument, IntField, EmbeddedDocumentField
from mongoengine.fields import StringField, ListField, DynamicField
from rdkit import Chem


class Identifiers(DynamicEmbeddedDocument):
    """IDs known for a molecule"""

    smiles = StringField(required=True)
    inchi = StringField(required=True)
    pubchem_id = IntField()


class MoleculeRecord(Document):
    """Defines whatever we know about a molecule"""

    # Identifiers
    key = StringField(min_length=27, max_length=27, required=True, primary_key=True, help_test='InChI key')
    identifier = EmbeddedDocumentField(Identifiers, help_text='Collection of identifiers which define the molecule')
    names = ListField(StringField(), help_text='Names this molecule is known by')

    @classmethod
    def from_identifier(cls, smiles: Optional[str] = None, inchi: Optional[str] = None):
        assert (smiles is not None) ^ (inchi is not None), "You must supply either smiles or inchi, and not both"

        # Load in the molecule
        if smiles is None:
            mol = Chem.MolFromInchi(inchi)
        else:
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f'Molecule failed to parse: {smiles if inchi is None else inchi}')

        # Create the object
        return cls(key=Chem.MolToInchiKey(mol), identifier=Identifiers(smiles=Chem.MolToSmiles(mol), inchi=Chem.MolToInchi(mol)))

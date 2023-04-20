from typing import Optional

from mongoengine import Document, DynamicEmbeddedDocument, EmbeddedDocument, IntField, EmbeddedDocumentField, MapField, FloatField
from mongoengine.fields import StringField, ListField
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

    # Characteristics
    subsets = ListField(StringField(), help_text='List of subsets this molecule is part of')
    property = MapField(FloatField(), help_text='Property values associated with this object')

    @classmethod
    def from_identifier(cls, smiles: Optional[str] = None, inchi: Optional[str] = None):
        if not ((smiles is None) ^ (inchi is None)):
            raise ValueError("You must supply either smiles or inchi, and not both")

        # Load in the molecule
        if smiles is None:
            mol = Chem.MolFromInchi(inchi)
        else:
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f'Molecule failed to parse: {smiles if inchi is None else inchi}')

        # Create the object
        return cls(key=Chem.MolToInchiKey(mol), identifier=Identifiers(smiles=Chem.MolToSmiles(mol), inchi=Chem.MolToInchi(mol)))


class Match(EmbeddedDocument):
    """Definition of a match"""

    key = StringField(help_text='InChI key of the matched molecule')
    name = StringField(help_text='Name of the molecule used in this match')


class Mention(Document):
    """Where a molecule was mentioned in text"""

    # We store one match per document
    filename = StringField(required=True, help_text='Name of the file containing the mention')
    line = IntField(required=True, unique_with='filename', help_text='Line number in the file')

    # Store the text for convenience
    text = StringField(required=True, help_text='Sentence containing the text')

    # The matches
    matches = ListField(EmbeddedDocumentField(Match), help_text='List of matches that appear in the dataset')

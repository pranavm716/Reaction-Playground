import pathlib
from enum import StrEnum, auto

import yaml
from rdkit import Chem
from rdkit.Chem import Mol


class MolClass(StrEnum):
    """
    Enum class representing the functional groups/substructures/classification of a molecule.
    """

    # Alcohols
    primary_alcohol = auto()
    secondary_alcohol = auto()
    tertiary_alcohol = auto()

    # Carbonyls
    aldehyde = auto()
    carboxylic_acid = auto()
    ketone = auto()
    ester = auto()

    # Nitrogen compounds
    primary_amine = auto()
    secondary_amine = auto()
    tertiary_amine = auto()

    primary_amide = auto()
    secondary_amide = auto()
    tertiary_amide = auto()

    nitrile = auto()

    # Alkyl bromides
    primary_alkyl_bromide = auto()
    secondary_alkyl_bromide = auto()
    tertiary_alkyl_bromide = auto()

    # Alkenes and alkynes
    alkene = auto()
    internal_alkyne = auto()
    terminal_alkyne = auto()

    # Other
    ether = auto()
    acid_chloride = auto()
    carbon_nucleophile = auto()

    def __new__(cls, *values):
        # Replace underscores with spaces for a more user-friendly representation
        value = str(*values)
        value_with_spaces = value.replace("_", " ")
        member = str.__new__(cls, value_with_spaces)
        member._value_ = value
        return member


# Mapping of the substructure name to its substructure patterns
SubstructureDict = dict[MolClass, list[Mol]]


def read_all_substructures_from_file(path: pathlib.Path) -> SubstructureDict:
    substructure_data: dict[str, list[str]] = yaml.safe_load(path.read_text())

    substructures: SubstructureDict = {}
    for name, smarts_list in substructure_data.items():
        name = MolClass(name)
        substructures[name] = [Chem.MolFromSmarts(smarts) for smarts in smarts_list]
    return substructures

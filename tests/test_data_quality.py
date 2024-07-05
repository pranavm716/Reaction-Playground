import io
from unittest.mock import patch

from rdkit import rdBase

from backend.config import ALL_REACTIONS_FILE_PATH, ALL_SUBSTRUCTURES_FILE_PATH
from backend.mol_classification import read_all_substructures_from_file
from backend.reaction import read_all_reactions_from_file


def test_data_is_configured_properly():
    """Ensures that rdkit does not throw any warnings when parsing all the reactions and the substructures."""
    rdBase.LogToPythonStderr()
    with patch("sys.stderr", new_callable=io.StringIO) as s:
        read_all_reactions_from_file(ALL_REACTIONS_FILE_PATH)
        read_all_substructures_from_file(ALL_SUBSTRUCTURES_FILE_PATH)

    assert not s.getvalue()

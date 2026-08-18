from pathlib import Path

import pytest

from tola.assembly.sanger_files import AssemblyYamlLocationError, find_yaml


@pytest.fixture
def files_dir():
    return Path(__file__).parent / "files"


@pytest.fixture
def latest_yaml(files_dir):
    return files_dir / (
        "Eutomostethus_luteiventris/assembly/draft/iyEutLute1.20250917/iyEutLute1.yaml"
    )


def test_find_yaml(files_dir, latest_yaml):
    working = files_dir / "Eutomostethus_luteiventris/working"
    yaml = find_yaml(working)
    assert yaml == latest_yaml

    with pytest.raises(AssemblyYamlLocationError):
        find_yaml(files_dir)


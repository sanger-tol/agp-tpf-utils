from pathlib import Path

import pytest

from tola.assembly.sanger_files import (
    AssemblyYaml,
    AssemblyYamlError,
    AssemblyYamlLocationError,
    find_yaml,
)


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

    # This might succeed if the repository is located in a directory tree
    # which has "assembly/branch" somewhere off any of its parent
    # directories!
    with pytest.raises(AssemblyYamlLocationError):
        find_yaml(files_dir)


def test_decon_file_names(latest_yaml):
    yml = AssemblyYaml(latest_yaml)
    assert yml.decontaminated_file_path(Path("base/mito.fa"), "mito") == Path(
        "base/mito.decontaminated.fa"
    )
    assert yml.decontaminated_file_path(Path("base/mito.fa.gz"), "mito") == Path(
        "base/mito.decontaminated.fa.gz"
    )
    assert yml.decontaminated_file_path(Path("base/mito.fasta.gz"), "mito") == Path(
        "base/mito.decontaminated.fa.gz"
    )
    with pytest.raises(AssemblyYamlError):
        yml.decontaminated_file_path(Path("base/mito.txt"), "mito")

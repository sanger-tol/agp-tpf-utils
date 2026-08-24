from pathlib import Path

import pytest

from tola.assembly.fragment import Fragment
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
    with pytest.raises(AssemblyYamlError):
        yml.decontaminated_file_path(Path("base/mito.fax.gz"), "mito")


def test_get_mito_assembly(latest_yaml):
    yml = AssemblyYaml(latest_yaml)

    mito = yml.mitochondrial_assembly
    assert mito is not None
    assert len(mito.scaffolds) == 1

    scffld = mito.scaffolds[0]
    assert scffld.name == "scaffold_MT_1"
    assert len(scffld.rows) == 1

    frag = scffld.rows[0]
    assert isinstance(frag, Fragment)
    assert frag.source == "mito"


def test_get_no_pltd_assembly(latest_yaml):
    yml = AssemblyYaml(latest_yaml)

    pltd = yml.chloroplast_assembly
    assert pltd is None




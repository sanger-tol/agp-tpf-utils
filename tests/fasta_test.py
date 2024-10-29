import pathlib

import pytest

from tola.fasta.index import index_fasta_bytes


def list_fasta_files():
    fasta_dir = pathlib.Path(__file__).parent / "fasta"
    for ff in fasta_dir.iterdir():
        if ff.suffix == '.fa':
            yield ff


@pytest.mark.parametrize("fasta_file", list_fasta_files())
def test_fai(fasta_file):
    fai_file = pathlib.Path(str(fasta_file) + ".fai")
    fai_str = fai_file.read_text()
    if not fai_file.exists():
        msg = f"Missing expected '.fai' file: {fai_file}"
        raise ValueError(msg)
    info = index_fasta_bytes(fasta_file)
    test_str = "".join(x.fai_row() for x in info)
    assert test_str == fai_str

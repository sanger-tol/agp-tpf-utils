import io
import pathlib

import pytest

# from tola.assembly.fragment import Fragment
# from tola.assembly.scaffold import Scaffold
from tola.fasta.index import FastaIndex, index_fasta_file, reverse_complement
from tola.fasta.stream import FastaStream


def list_fasta_files():
    fasta_dir = pathlib.Path(__file__).parent / "fasta"
    for ff in fasta_dir.iterdir():
        if ff.suffix == ".fa":
            yield ff


@pytest.mark.parametrize("fasta_file", list_fasta_files())
def test_fai(fasta_file):
    idx = FastaIndex(fasta_file)
    idx.load_index()
    idx_dict, asm = index_fasta_file(fasta_file)
    assert idx_dict == idx.index
    idx.load_assembly()
    asm.header = idx.assembly.header = []
    assert str(asm) == str(idx.assembly)


def test_stream_fetch():
    fasta_file = pathlib.Path(__file__).parent / "fasta/test.fa"
    fai = FastaIndex(fasta_file, buffer_size=7)
    fai.load_index()
    fai.load_assembly()
    out = io.BytesIO()
    fst = FastaStream(out, fai, gap_character=b"n")
    fst.write_assembly(fai.assembly)

    print(out.getvalue().decode(), end="")

    ref_bytes = fasta_file.read_bytes().replace(b"\r", b"").replace(b"\n", b"")
    fst_bytes = out.getvalue().replace(b"\n", b"")

    assert len(ref_bytes) == len(fst_bytes)
    assert ref_bytes == fst_bytes


def test_revcomp():
    seq = b"ACGTRYMKSWHBVDNacgtrymkswhbvdn"
    assert reverse_complement(seq) == b"nhbvdwsmkryacgtNHBVDWSMKRYACGT"

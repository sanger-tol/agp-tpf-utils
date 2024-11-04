import io
import pathlib

import pytest

# from tola.assembly.fragment import Fragment
# from tola.assembly.scaffold import Scaffold
from tola.fasta.index import FastaIndex, index_fasta_file
from tola.fasta.simple import FastaSeq, reverse_complement
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


def test_simple_fasta_bytes():
    name = "test"
    desc = "A test sequence"
    seq = b"n" * 60
    ref_str = f">{name} {desc}\n{seq.decode()}\n"

    fst = FastaSeq(name, seq, desc)
    assert str(fst) == ref_str
    assert fst.fasta_bytes() == ref_str.encode()


@pytest.mark.parametrize("buf_size", [5, 7, 100, 200])
def test_stream_fetch(buf_size):
    fasta_file = pathlib.Path(__file__).parent / "fasta/test.fa"
    ref_fai = FastaIndex(fasta_file)
    ref_fai.load_index()

    # Check we have the first and last sequence
    assert ref_fai.index.get('RAND-001')
    assert ref_fai.index.get('RAND-100')

    ref_io = io.BytesIO()
    for seq in ref_fai.all_fasta_seq():
        ref_io.write(seq.fasta_bytes())
    ref_bytes = ref_io.getvalue()

    fai = FastaIndex(fasta_file, buffer_size=buf_size)
    fai.load_index()
    fai.load_assembly()

    out = io.BytesIO()
    fst = FastaStream(out, fai, gap_character=b"n")
    fst.write_assembly(fai.assembly)
    fst_bytes = out.getvalue()

    # Decode bytes to string so that pytest diff works
    assert ref_bytes.decode() == fst_bytes.decode()


def test_revcomp():
    seq = b"ACGTRYMKSWHBVDNacgtrymkswhbvdn"
    assert reverse_complement(seq) == b"nhbvdwsmkryacgtNHBVDWSMKRYACGT"

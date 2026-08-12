import io
import pathlib

import pytest

from tola.assembly.assembly import Assembly
from tola.fasta.index import FastaIndex, index_fasta_file
from tola.fasta.simple import FastaSeq, reverse_complement
from tola.fasta.stream import FastaCollection, FastaStream


def list_fasta_files():
    fasta_dir = pathlib.Path(__file__).parent / "fasta"
    return [ff for ff in fasta_dir.iterdir() if ff.suffix == ".fa"]


@pytest.mark.parametrize("fasta_file", list_fasta_files())
def test_fai(fasta_file):
    idx = FastaIndex(fasta_file)
    idx.load_index()
    idx_dict, asm = index_fasta_file(fasta_file)
    assert idx_dict == idx.index
    idx.load_assembly()
    asm.header = idx.assembly.header = []
    assert str(asm) == str(idx.assembly)


def test_revcomp():
    seq = b"ACGTRYMKSWHBVDNacgtrymkswhbvdn"
    assert reverse_complement(seq) == b"nhbvdwsmkryacgtNHBVDWSMKRYACGT"


def test_simple_fasta_bytes():
    name = "test"
    desc = "A test sequence"
    seq = b"atcg" * 15
    rev = b"cgat" * 15
    ref_str = f">{name} {desc}\n{seq.decode()}\n"
    rev_ref_str = f">{name} {desc}\n{rev.decode()}\n"

    fst = FastaSeq(name, seq, desc)
    assert str(fst) == ref_str
    assert ref_str.encode() == fst.fasta_bytes()

    rev_fst = fst.rev_comp()
    assert str(rev_fst) == rev_ref_str
    assert rev_ref_str.encode() == rev_fst.fasta_bytes()


@pytest.mark.parametrize("buf_size", [5, 7, 100, 200])
def test_stream_fetch(buf_size):
    """
    Tests with small buffer sizes so that the chunking code is exercised.
    """

    fasta_file = pathlib.Path(__file__).parent / "fasta/test.fa"
    ref_fai = FastaIndex(fasta_file)
    ref_fai.load_index()

    # Check we have the first and last sequence
    assert ref_fai.index.get("RAND-001")
    assert ref_fai.index.get("RAND-100")

    ref_io = io.BytesIO()
    for seq in ref_fai.all_fasta_seq():
        ref_io.write(seq.fasta_bytes())
    ref_bytes = ref_io.getvalue()

    fai = FastaIndex(fasta_file, buffer_size=buf_size)
    fai.load_index()
    fai.load_assembly()
    coll = FastaCollection(fai)

    out = io.BytesIO()
    fst = FastaStream(out, coll, gap_character=b"n")
    fst.write_assembly(fai.assembly)
    fst_bytes = out.getvalue()

    # Decode bytes to string so that pytest diff works
    assert ref_bytes.decode() == fst_bytes.decode()

    # Make reverse-complement reference
    rev_ref_io = io.BytesIO()
    for seq in ref_fai.all_fasta_seq():
        rev_ref_io.write(seq.rev_comp().fasta_bytes())
    rev_ref_bytes = rev_ref_io.getvalue()

    # Test reverse-complement of assembly
    rev_out = io.BytesIO()
    rev_fst = FastaStream(rev_out, coll, gap_character=b"n")
    for scffld in fai.assembly.scaffolds:
        scffld_rev = scffld.reverse()
        rev_fst.write_scaffold(scffld_rev)
    rev_fst_bytes = rev_out.getvalue()
    assert rev_ref_bytes.decode() == rev_fst_bytes.decode()


def test_multi_file_collection():
    fasta_dir = pathlib.Path(__file__).parent / "fasta"

    fa1 = FastaIndex(fasta_dir / "test.fa")
    fa1.load_index()
    fa1.load_assembly()

    fa2 = FastaIndex(fasta_dir / "test_other.fa")
    fa2.load_index()
    fa2.load_assembly()

    coll = FastaCollection()
    coll.add_faidx(fa1)
    coll.add_faidx(fa2, "othr")

    ref_io = io.BytesIO()
    for seq in (
        fa1.get_fasta_seq("RAND-011"),
        fa2.get_fasta_seq("RAND-003"),  # From second index
        fa1.get_fasta_seq("RAND-090"),
    ):
        ref_io.write(seq.fasta_bytes())
    ref_str = ref_io.getvalue().decode()

    new_asm = Assembly(name="test-mixed")
    new_asm.add_scaffold(fa1.assembly.scaffolds[10])
    new_asm.add_scaffold(fa2.assembly.scaffolds[2])
    new_asm.add_scaffold(fa1.assembly.scaffolds[89])

    new_out = io.BytesIO()
    fst = FastaStream(new_out, coll, gap_character=b"n")
    fst.write_assembly(new_asm)
    new_str = new_out.getvalue().decode()
    assert new_str == ref_str



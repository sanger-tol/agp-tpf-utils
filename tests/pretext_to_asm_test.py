import difflib
import io
import re
import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from tola.assembly.scripts.pretext_to_asm import cli, parse_output_file


def test_parse_output_file():
    assert parse_output_file(Path("xx/aaa.2.fa")) == (
        "FASTA",
        Path("xx"),
        "aaa",
        "2",
        ".fa",
    )
    assert parse_output_file(Path("xx/aaa.fasta")) == (
        "FASTA",
        Path("xx"),
        "aaa",
        "1",
        ".fasta",
    )
    assert parse_output_file(Path("xx/aaa.agp")) == (
        "AGP",
        Path("xx"),
        "aaa",
        "1",
        ".agp",
    )


def list_example_assemblies():
    data_dir = Path(__file__).parent / "data"
    for xd in data_dir.iterdir():
        if xd.is_dir():
            yield xd


@pytest.mark.parametrize("specimen_dir", list_example_assemblies())
def test_assembly(specimen_dir):
    specimen = specimen_dir.name
    version = ""
    if m := re.search(r"_(\d+)$", specimen):
        version = "." + m.group(1)
        specimen = specimen[: -len(version)]
    input_tpf = f"{specimen}-input{version}.tpf"
    pretext_agp = f"{specimen}-pretext{version}.agp"
    assert (specimen_dir / input_tpf).exists()
    assert (specimen_dir / pretext_agp).exists()

    output_tpf = f"{specimen}-pretext-to-tpf{version}.tpf"
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        args = (
            "--assembly",
            specimen_dir / input_tpf,
            "--pretext",
            specimen_dir / pretext_agp,
            "--output",
            tmp_path / output_tpf,
            "--write-log",
        )
        runner = CliRunner()
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        for spec_path in specimen_dir.iterdir():
            spec_file = spec_path.name
            if spec_file in (input_tpf, pretext_agp) or spec_file.startswith("."):
                # Skip the input and hidden files
                continue
            test_path = tmp_path / spec_file
            print(f"Specimen dir file: {spec_file}")
            assert test_path.exists()
            diff_files(specimen, spec_path, test_path)


def diff_files(specimen, spec_path, test_path):
    ctx_diff = difflib.context_diff(
        spec_path.open("r").readlines(),
        test_path.open("r").readlines(),
        fromfile=f"{specimen}/{spec_path.name}",
        tofile=f"test/{test_path.name}",
    )
    str_io = io.StringIO()
    str_io.writelines(ctx_diff)
    ctx_str = str_io.getvalue()
    if ctx_str:
        print(ctx_str)
    assert ctx_str == ""

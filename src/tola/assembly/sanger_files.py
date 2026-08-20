"""
For finding files in the Wellcome Sanger Institute's genome assembly filesystem.
"""

import logging
import sys
from functools import cached_property
from pathlib import Path

import yaml

from tola.assembly.assembly import Assembly
from tola.fasta.index import FastaIndex
from tola.fasta.stream import FastaCollection

log = logging.getLogger(__name__)


class AssemblyYamlLocationError(Exception):
    """The YAML file for the draft assembly was not found"""


def find_yaml(search_dir: Path) -> Path:
    abs_dir = search_dir.absolute().resolve()

    asm_draft = None
    for d in abs_dir.parents:
        maybe_asm = d / "assembly/draft"
        if maybe_asm.exists():
            asm_draft = maybe_asm
            break

    if not asm_draft:
        msg = (
            "Failed to find an 'assembly/draft' directory"
            f" in any parent directory of '{abs_dir}'"
        )
        raise AssemblyYamlLocationError(msg)

    yaml_files = sorted(asm_draft.glob("*/*.yaml"))
    if not yaml_files:
        msg = f"Failed to find a YAML file in any sub directory of '{asm_draft}'"
        raise AssemblyYamlLocationError(msg)

    yaml = yaml_files.pop()
    if yaml_files:
        log.warning(
            "\n  ".join(
                (
                    f"Found multiple YAML files. Chose '{yaml}' and ignored:",
                    *[f"'{x}'" for x in yaml_files],
                ),
            )
        )

    return yaml


class AssemblyYamlError(Exception):
    """Error in expectied contents of genome assembly YAML file"""


class AssemblyYaml:
    def __init__(self, yaml_file: Path):
        self.__yaml_dir = yaml_file.parent
        self.__load_yaml(yaml_file)
        self.__fasta_index_list = []

    def __load_yaml(self, yaml_file: Path) -> None:
        yaml_dict = yaml.safe_load(yaml_file)
        if not isinstance(yaml_dict, dict):
            msg = f"Expected a dict from '{yaml_file}' but got {yaml_dict!r}"
            raise AssemblyYamlError(msg)
        self.__yaml: dict[str, str] = yaml_dict

    def add_indexes_to_collection(self, coll: FastaCollection):
        for fai in self.__fasta_index_list:
            coll.add_faidx(fai)

    @cached_property
    def mitochondrial_assembly(self):
        return self.index_fasta("mito")

    @cached_property
    def chloroplast_assembly(self):
        return self.index_fasta("plastid")

    @cached_property
    def haplotigs_assembly(self):
        return self.index_fasta("haplotigs")

    def index_fasta(self, name: str) -> Assembly | None:
        file = self.__yaml.get(name)
        if file is None:
            return None

        file_path = Path(file)
        if not file_path.is_absolute:
            file_path = self.__yaml_dir / file_path
        if not file_path.exists():
            msg = "No such FASTA file for {name!r}: '{file_path}'"
            raise AssemblyYamlError(msg)
        fai = FastaIndex(file_path, source=name)
        self.__fasta_index_list.append(fai)
        return fai.assembly


if __name__ == "__main__":
    for working_dir in sys.argv[1:]:
        print(find_yaml(Path(working_dir)))

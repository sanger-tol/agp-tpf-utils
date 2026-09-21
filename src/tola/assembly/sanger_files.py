"""
For finding files in the Wellcome Sanger Institute's genome assembly filesystem.
"""

import logging
import re
import sys
from functools import cached_property
from pathlib import Path
from typing import Any

import yaml

from tola.assembly.assembly import Assembly
from tola.fasta.index import FastaIndex
from tola.fasta.stream import FastaCollection

log = logging.getLogger(__name__)


class AssemblyYamlLocationError(Exception):
    """The YAML file for the draft assembly was not found"""


def find_yaml(search_dir: Path, branch_dir: str | Path = "assembly/draft") -> Path:
    """
    Looks in each parent directory of `search_dir` for the `branch_dir`
    structure, and then looks for any `*.yaml` files in each sub-directory of
    the first `branch_dir` found. *i.e.*

    ```python
    search_dir = "base/working"
    branch_dir = "assembly/draft"
    yaml_files = "base/assembly/draft/*/*.yaml"
    ```

    When multiple YAML files are found, prints a warning and returns the one
    that sorts lexically last.

    The returned file is an absolute `Path`.  Throws a
    `AssemblyYamlLocationError` exception on failure.
    """
    abs_dir = search_dir.absolute().resolve()

    asm_draft = None
    for d in abs_dir.parents:
        maybe_asm = d / branch_dir
        if maybe_asm.exists():
            asm_draft = maybe_asm
            break

    if not asm_draft:
        msg = (
            f"Failed to find an '{branch_dir}' directory"
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
                    "Found multiple YAML files. Chose:",
                    f"'{yaml}'\nIgnored:",
                    *[f"'{x}'" for x in yaml_files],
                ),
            )
        )

    return yaml


class AssemblyYamlError(Exception):
    """Error in expectied contents of genome assembly YAML file"""


class AssemblyYaml:
    def __init__(self, yaml_file: Path):
        self.__file = yaml_file
        self.__yaml_dir = yaml_file.parent.absolute()
        self.__load_yaml(yaml_file)
        self.__fasta_index_list = []

    @property
    def file(self):
        return self.__file

    def __load_yaml(self, yaml_file: Path) -> None:
        yaml_dict = yaml.safe_load(yaml_file.open())
        if not isinstance(yaml_dict, dict):
            msg = f"Expected a dict from '{yaml_file}' but got {yaml_dict!r}"
            raise AssemblyYamlError(msg)
        self.__yaml: dict[str, Any] = yaml_dict

    def get_path(self, key: str) -> Path:
        val = self.__yaml.get(key)
        if val is None:
            msg = f"No such value {key!r}"
            raise AssemblyYamlError(msg)
        if not isinstance(val, str):
            msg = f"Expecting string for {key!r} but got: {val!r}"
            raise AssemblyYamlError(msg)

        path = Path(val)
        return path if path.is_absolute() else self.__yaml_dir.joinpath(path).resolve()

    def add_indexes_to_collection(self, coll: FastaCollection):
        for fai in self.__fasta_index_list:
            coll.add_faidx(fai)

    @cached_property
    def mitochondrial_assembly(self) -> Assembly | None:
        return self.__get_asm_and_name_scaffolds("mito", "MT", rank=3, localised=True)

    @cached_property
    def chloroplast_assembly(self) -> Assembly | None:
        return self.__get_asm_and_name_scaffolds("plastid", "Pltd", rank=3, localised=True)

    def __get_asm_and_name_scaffolds(
        self,
        yaml_key,
        prefix,
        *,
        rank: int = 4,
        localised: bool = False,
    ) -> Assembly | None:
        asm = self.index_fasta(yaml_key)
        if not asm:
            return None

        ### Mark linear organelle genomes here? ###

        multi_flag = len(asm.scaffolds) > 1

        # Give each scaffold a name and set its rank
        for i, scffld in enumerate(asm.scaffolds, start=1):
            scffld.name = f"scaffold_{prefix}_{i}"
            scffld.chr_name = f"{prefix}-{i}" if multi_flag else prefix
            scffld.rank = rank
            scffld.localised = localised
        return asm

    @cached_property
    def haplotigs_assembly(self) -> Assembly | None:
        return self.index_fasta("haplotigs")

    def index_fasta(self, name: str) -> Assembly | None:
        file_name = self.__yaml.get(name)
        if file_name is None:
            return None

        file_path = Path(file_name)
        if not file_path.is_absolute():
            file_path = self.__yaml_dir / file_path

        decon_path = self.decontaminated_file_path(file_path, name)

        # Do we have the "decontaminated" version of the file?
        if not decon_path.exists():
            msg = f"No such FASTA file for {name!r}: '{decon_path}'"
            raise AssemblyYamlError(msg)

        # Index the FASTA with a source namespace
        fai = FastaIndex(decon_path, source=name)
        fai.run_indexing()  # Only run indexing.  Don't write files
        self.__fasta_index_list.append(fai)
        return fai.assembly

    def decontaminated_file_path(self, file_path: Path, name: str) -> Path:
        """
        Construct *e.g.* `draft/mito.decontaminated.fa.gz` from
        `draft/mito.fa.gz`.  Input file can end with `.fasta.gz` but the
        decontaminated file will always end `.fa.gz`.  If the input
        `file_path` has no `.gz` suffix, the decontaminated path returned
        will not have it either.
        """
        base_gz = re.split(r"\.fa(?:sta)?\b", file_path.name)
        if len(base_gz) != 2:
            msg = f"Unexpected FASTA file name format for {name!r}: '{file_path}'"
            raise AssemblyYamlError(msg)
        base, gz = base_gz
        return file_path.parent / f"{base}.decontaminated.fa{gz.lower()}"


if __name__ == "__main__":
    for working_dir in sys.argv[1:]:
        print(find_yaml(Path(working_dir)))

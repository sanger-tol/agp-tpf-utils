"""
For finding files in the Wellcome Sanger Institute's genome assembly filesystem.
"""

import logging
import sys
from pathlib import Path

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
                f"Found multiple YAML files. Chose '{yaml}' and ignored:",
                *[f"'{x}'" for x in yaml_files],
            )
        )

    return yaml


if __name__ == "__main__":
    for working_dir in sys.argv[1:]:
        print(find_yaml(Path(working_dir)))

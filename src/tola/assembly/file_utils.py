import logging
import sys
from pathlib import Path
from typing import IO, Any

import click
from zlib_ng import gzip_ng_threaded

log = logging.getLogger(__name__)


def get_output_filehandle(path: Path, clobber: bool, mode: str = "") -> IO[Any]:
    op = "Overwrote" if path.exists() else "Created"

    # Must choose binary output mode if output is gzip compressed
    gz = path.suffix == ".gz"
    if gz:
        mode = "b"
    mode = "w" + mode if clobber else "x" + mode

    try:
        out_fh = (
            gzip_ng_threaded.open(path, mode, compresslevel=6, threads=2)
            if gz
            else path.open(mode)
        )
    except FileExistsError:
        log.error(f"Output file '{path}' already exists")
        sys.exit(1)
    click.echo(f"{op:>11}: '{path}'", err=True)
    return out_fh

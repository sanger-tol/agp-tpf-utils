from pathlib import Path

import click

from tola.fasta.index import FastaIndex


@click.command
@click.argument(
    "fasta_files",
    nargs=-1,
    required=True,
    type=click.Path(
        path_type=Path,
        exists=True,
        readable=True,
        dir_okay=False,
    ),
    metavar="FASTA_FILE...",
)
def cli(fasta_files):
    """
    Index all the FASTA files provided on the command line. An ".agp"
    and ".fai" file will be created (or recreated if existing index files are
    older than the FASTA file).

    Each FASTA file can be gzip compressed, in which case they are expected to
    end in the suffix ".gz", and a ".gzidx" will be created too.
    """

    for ff in fasta_files:
        idx = FastaIndex(ff)
        if not idx.check_for_index_files():
            idx.run_indexing()
            idx.write_index()
            idx.write_assembly()
            idx.write_gzip_index()



# Utilities for Tree of Life AGP and TPF Assembly Files

Code for working with AGP and TPF files as used within the Tree of Life
project, where the combination of long read sequencing and Hi-C data is used
to produce whole genome assemblies. It is not therefore intended to cover the
full range of AGP and TPF syntax.

## Scripts

Added to your `PATH` if the suggested development venv is set up. Run with
`--help` for usage.

### [`asm-format`](src/tola/assembly/scripts/asm_format.py)

Parses and reformats AGP and TPF files, converting into either format.

### [`pretext-to-tpf`](src/tola/assembly/scripts/pretext_to_tpf.py)

Takes the AGP file output by
[PretextView](https://github.com/wtsi-hpag/PretextView)
and creates TPF files containing precise coordinates of the curated assembly.

## File Formats

Both TPF and AGP file formats described here contain the same information. AGP
is the more appropriate format to use, since it was designed for sequence
assembly coordinates, whereas TPF was for listing (cosmid, fosmid, YAC or
BAC) clones and their accessions in the order that they were tiled to build a
chromosome.

### AGP

Each line in the
[AGP v2.1 specification](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/)
contains 9 tab delimited columns. Of these columns:

- **DNA Sequence**
    - **column 5** the "component_type" contains `W` in our assemblies,
        meaning a contig from Whole Genome Shotgun (WGS) sequencing.
    - **columns 10 and greater** are extra tag metadata columns not included
        in the AGP v2.1 specification. (See below for their possible
        values.)
- **Gaps**
    - **column 5** the "component_type" contains `U` in our assemblies, for a
        gap of unknown length. (The other gap type `N` is for gaps of known
        length.)
    - **column 6** The default length in the specification for `U` gaps is 100
        base pairs, but we use 200 bp gaps, as produced by
        [yahs](https://github.com/sanger-tol/yahs)
    - **column 7** has `scaffold`, signifying a gap between two contigs in a
        scaffold.
    - **column 8** has `yes`, signifying that there is evidence of linkage
        between the sequence data on either side of the gap.

#### Tags

Single words appended in tab-delimted columns beyond column 9, they can
contain:

- `Contaminant`
- `Haplotig` for haplotype-specific contigs.
- Haplotypes:
  - `Hap1`, `Hap2`…
- `Painted` where fragment has Hi-C contacts.
- `Target`
- `Unloc` are fragments attached to chromosomes but unlocalised within them.
- Sex Chromosomes:
  - `U`
  - `V`
  - `W` or `W1`, `W2`…
  - `X` or `X1`, `X2`…
  - `Y` or `Y1`, `Y2`…
  - `Z` or `Z1`, `Z2`…
- [B Chromosomes](https://en.wikipedia.org/wiki/B_chromosome):
  - `B1`, `B2`, `B3`…

### TPF

Our TPF files are highly diverged from the
[original specification](https://www.ncbi.nlm.nih.gov/projects/genome/assembly/TPF_Specification_v1.4_20110215.pdf).

- We incorporate assembly coordinates, which was not the purpose of TPF files.
- We do not necessarily include any `##` header lines, which were mandatory in
  the original specification.
- **DNA Sequence**
    - **column 1** the "accession" is always `?` since the components of our
        assemblies are not accessioned.
    - **column 2** the "clone name" does not contain a clone name, but
        contains the name of scaffold fragment or whole scaffold, with the
        format: `<name>:<start>-<end>` *i.e.* assembly coordinates.
    - **column 3** the "local contig identifier" now contains the name of the
        scaffold each sequence fragment belongs to. Each TPF file used to
        contain a single chromosome, but we put a whole genome into a single
        file, and this column groups the fragments into chromosomes /
        scaffolds.
    - **column 4** which in the original specification was used for
        indicating `CONTAINED` or `CONTAINED_TURNOUT` clones now holds
        assembly strand information, either `PLUS` or `MINUS`.
- **Gaps**
    - **column 2** is `TYPE-2`, which meant a gap between two clones
    - **column 3** length, using our default of 200 bp.

## Development Setup

In your cloned copy of the git repository:

```sh
python3 -m venv --prompt asm-utils venv
source venv/bin/activate
pip install --upgrade pip
pip install --editable .
```

An alias such as this:

```sh
alias atu="cd $HOME/git/agp-tpf-utils && source ./venv/bin/activate"
```

in your shell's `.*rc` file (*e.g.* `~/.bashrc` for `bash` or `~/.zshrc` for
`zsh`) can be convenient.

### Reinstalling Development Environment

Some changes, such as adding a new command line script to
[`pyproject.toml`](pyproject.toml), require the development environment to be
reinstalled:

```sh
pip uninstall tola-agp-tpf-utils
pip install --editable .
hash -r
```

## Running Tests

Tests, located in the [`tests/`](tests) directory, are run with the `pytest`
command from the project root.

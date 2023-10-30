
# Assembly Utilities for Tree of Life

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
cd $HOME/git/assembly-utils && source ./venv/bin/activate
```

in your shell's `.*rc` file (*e.g.* `~/.bashrc` for `bash` or `~/.zshrc` for
`zsh`) can be convenient.

### Reinstalling Development Environment

Some changes, such as adding a new command line tool to
[`pyproject.toml`](pyproject.toml), require the development environment to be
reinstalled:

```sh
pip uninstall tola-assembly-utils
pip install --editable .
```

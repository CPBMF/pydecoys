[![Issues](https://img.shields.io/badge/Issues-critical?logo=github)](https://github.com/CPBMF/pydecoys/issues)
[![CI](https://github.com/CPBMF/pydecoys/actions/workflows/tests.yaml/badge.svg)](https://github.com/CPBMF/pydecoys/actions/workflows/tests.yaml)
[![documentation](https://github.com/CPBMF/pydecoys/actions/workflows/documentation.yaml/badge.svg)](https://github.com/CPBMF/pydecoys/actions/workflows/documentation.yaml)
[![Coverage](https://CPBMF.github.io/pydecoys/_static/coverage.svg)](https://pypi.org/project/pytest-cov/)
![Image](https://img.shields.io/badge/Python-3.12+-FFD43B?logo=python&logoColor=blue)
[![Image](https://img.shields.io/badge/License-GPLv3+-red?logo=gplv3)](https://github.com/CPBMF/pydecoys/blob/main/LICENSE)

# PyDecoys

Welcome to the GitHub page of PyDecoys:
Proteomics decoy utilities for Python!

Check the full documentation here: <https://cpbmf.github.io/pydecoys/>.

## About

PyDecoys is a bioinformatics Python package and application for generating decoy
proteins from target proteins.

Decoy generation is an essential step of Proteomics workflows: decoy proteins are
included in the search space as necessarily false hits; assuming decoy matches are as
likely as target false-positive matches, decoy/target match ratios allow for the
estimation of false-positive error ratios.

PyDecoys aims to facilitate this decoy generation step in a way that's easily
integratable to existing Proteomics workflows. It provides:

1. A CLI app to easily generate a decoy fasta file from a target fasta file
1. A comprehensive Python API with IO, Generator and single data functions
1. Easy implementation of custom decoy strategies
1. Full integration with [Biopython]

### Decoy strategies available

The following are the implemented strategies:

- Reverse proteins
- Shuffle proteins
- Randomize proteins
- Rerverse enzymatic peptides
- Shuffle enzymatic peptides
- Randomize enzymatic peptides

The most used proteases are available:

- Trypsin
- Strict trypsin
- Arg-C
- Asp-N
- Glu-C
- Lys-C
- Lys-N
- Chymotrypsin
- Pepsin-A
- CNBr

A full tutorial on how to set new proteases and decoy strategies can be found
at the [documentation].

## Getting started

PyDecoys is available at PyPI and can be easily set up.

### Dependencies

PyDecoys requires Python 3.12 or superior.

### Installation

If you only care about the CLI, PyDecoys can be easily installed or uninstalled
via `pipx`:

```sh
# Install current version
pipx install git+https://github.com/CPBMF/pydecoys@v0.2.1

# Uninstall
pipx uninstall pydecoys
```

`pipx` should automatically install PyDecoys and make it globally available on
PATH. You can confirm it's available by running ``pydecoys -v``. In case it
isn't, run:

```sh
pipx ensurepath
```

This will ensure all `pipx` scripts are available. The `pipx` app itself can be
installed via `pip`:

```sh
pip install --user pipx
```

In case `pipx` didn't automatically install [Biopython], run:

```sh
pipx inject pydecoys biopython
```

For API usage, you can use `pip` or other package managers:

```sh
# Install
pip install git+https://github.com/CPBMF/pydecoys@v0.2.1

# Uninstall
pip uninstall pydecoys
```

#### Installing without Biopython

To install PyDecoys without [Biopython], use the ``--no-deps`` flag:

```sh
pip install git+https://github.com/CPBMF/pydecoys@v0.2.1 --no-deps
```

If you change your mind later, simply run ``pip install biopython``.
[Biopython]'s dependency group is `biopython`.

Note that IO functions and the CLI app **aren't available without Biopython**.

## Usage

PyDecoys can be used both as a CLI app and through its Python interface.

### CLI

To generate a reversed-proteins decoy from a given `targets.fasta`, run:

```sh
pydecoys targets.fasta -o rev.fasta
```

You can change the decoy strategy with the `-s/--strategy` flag:

```sh
pydecoys targets.fasta -o shuf.fasta -s shuffle
```

The `input` defaults to `stdin`, so you can redirect the output of another
command to `pydecoys`.

```sh
cat database_1.fasta database_2.fasta | pydecoys -o rev_1_2.fasta
```

The same can be done with the output:

```sh
pydecoys targets.fasta | grep some_seq
```

Run with `-h/--help` to get a help message. Run with `-v/--version` to see the
version. A full list of CLI flags can be found at the [documentation].

### API

PyDecoys has a comprehensive API for integration in Python proteomics
workflows.

```py
import pydecoys
```

To get decoy sequences from a target fasta file, do:

```py
decoys = pydecoys.from_fasta('path/to/targets.fasta', strategy='reverse')
```
You can also write a target dataset to a fasta file:

```py
from Bio import SeqIO

# `to_fasta` uses the Biopython interface
targets = SeqIO.read('path/to/targets.fasta', format='fasta')
pydecoys.to_fasta(targets, 'path/to/decoys.fasta', strategy='reverse')
```

For convenience, you can directly pass a path or file handle to `to_fasta`:

```py
pydecoys.to_fasta('path/to/targets.fasta', 'path/to/decoys.fasta', strategy='reverse')
```

Alongside IO functions, there also are some generator functions for decoy
generation. To get decoys from [Biopython] `SeqRecord` objects:

```py
targets: list[SeqRecord]
decoys = pydecoys.from_SeqRecords(targets, strategy='reverse')
```

You can also opt-out of [Biopython] and either use a `tuple[str, str]` iterable
where the first item is the seqid and the second item is the aa sequence, or
just a `str` iterable representing the aa sequence:

```py
# With tuples...
decoys = pydecoys.from_tuples(tuple_targets, strategy='reverse')

# With strings...
decoys = pydecoys.from_seqs(str_targets, strategy='reverse')
```

Both functions **can** use [Biopython] `Seq` or `MutableSeq` objects as aa
sequence representations, but they don't need [Biopython] as a dependency.

Finally, creating a new strategy can be as simple as writing a function:

```py
def custom_decoys(target):
    # Creating decoy...
    return decoy

decoys = pydecoys.from_fasta('path/to/targets.fasta', custom_decoys)

# You can also register a `str` key for it:
pydecoys.register('custom-decoys', custom_decoys)
decoys = pydecoys.from_fasta('path/to/targets.fasta', 'custom-decoys')
```

For more options and examples, as well as a guide on the `strategies` module,
refer to the [documentation].

## Reporting issues

You can report issues through GitHub: <https://github.com/CPBMF/pydecoys/issues>.

Provide full context of how the issue happened, as well as how to reproduce
it. If possible, please provide a code snippet.

## Authorship and affiliations

This software was developed by Bruno Maestri A Becker as part of a research
project under the supervision of Professor Cristiano Valim Bizarro, PhD, at
Centro de Pesquisas em Biologia Molecular e Funcional (CPBMF) and Instituto
Nacional de Ciência e Tecnologia em Tuberculose (INCT-TB), Pontifícia
Universidade Católica do Rio Grande do Sul (PUCRS), Porto Alegre, Brazil.

ORCID™ iDs:
- [![ORCID][iD]](https://orcid.org/0009-0002-2338-5997) `Bruno Maestri A Becker`
- [![ORCID][iD]](https://orcid.org/0000-0002-2609-8996) `Cristiano Valim Bizarro`

[Biopython]: https://github.com/biopython/biopython
[documentation]: https://cpbmf.github.io/pydecoys/
[iD]: https://img.shields.io/badge/ORCID-iD-%23A6CE39?logo=orcid

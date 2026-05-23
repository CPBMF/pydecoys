# Copyright (C) 2026 CPBMF and INCT-TB, PUCRS, Porto Alegre, Brazil
# Copyright (C) 2026 Bruno Maestri A Becker
#
# This file is part of PyDecoys.
#
# PyDecoys is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PyDecoys is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PyDecoys. If not, see <https://www.gnu.org/licenses/>.

"""
Provides
    1. A CLI app to easily generate a decoy fasta file from a target fasta file
    1. A varied set of decoy generation strategies
    1. A comprehensive Python API
    1. Easy implementation of custom decoy strategies
    1. Full integration with Biopython

Available functionality
-----------------------
:func:`from_fasta`
    Return a decoy generator from a target fasta (requires `Biopython`)
:func:`to_fasta`
    Write a decoy fasta from a target protein set or fasta (requires `Biopython`)
:func:`from_SeqRecords`
    Return a decoy generator from a target set of `SeqRecord` (requires `Biopython`)
:func:`from_seqs`
    Return a decoy generator from a target set of :obj:`SeqLike`
:func:`from_tuples`
    Return a decoy generator from a target set of `tuple`
:func:`SeqRecord_as_decoy`
    Return a decoy `SeqRecord` from a given `SeqRecord` (requires `Biopython`)
:func:`seq_as_decoy`
    Return a decoy :obj:`SeqLike` from a given :obj:`SeqLike`
:func:`tuple_as_decoy`
    Return a decoy `tuple` from a given `tuple`
:func:`register`
    Register a custom decoy strategy on the API
:obj:`SeqLike`
    Custom type ``(str | Seq | MutableSeq)``
:mod:`strategies`
    API for custom decoy generation strategies

Avaliable decoy strategies
--------------------------
Each decoy strategy is specified by a lowercase string. The
:class:`strategies.DecoyGenerator` fns used for each strategy are available
under :mod:`strategies`. New strategies can be created by following the
:class:`strategies.DecoyGenerator` protocol or by instantiating
:class:`strategies.PseudoReverseRule` or :class:`strategies.PseudoShuffleRule`.

- reverse:                   Reverse protein
- reverse-keepn:             Reverse protein, except N-terminal aa
- reverse-keepc:             Reverse protein, except C-terminal aa
- reverse-keepterm:          Reserse protein, except terminal aas
- shuffle:                   Shuffle protein
- shuffle-keepn:             Shuffle protein, except N-terminal aa
- shuffle-keepc:             Shuffle protein, except C-terminal aa
- shuffle-keepterm:          Shuffle protein, except terminal aas
- reversepep-<enzyme>:       Pseudo-reverse enzymatic fragments
- reversepep-<enzyme>-keepn: Pseudo-reverse enzymatic fragments, except N-terminal aa
- shufflepep-<enzyme>:       Pseudo-shuffle enzymatic fragments
- shufflepep-<enzyme>-keepn: Pseudo-shuffle enzymatic fragments, except N-terminal aa

You can check the full list of proteases and how to add more at the documentation.
"""  # noqa: W505, E501

from importlib.metadata import version


__version__ = version('pydecoys')

__all__ = [
    'from_fasta',
    'to_fasta',
    'from_SeqRecords',
    'from_seqs',
    'from_tuples',
    'SeqRecord_as_decoy',
    'seq_as_decoy',
    'tuple_as_decoy',
    'register',
    'strategies',
    'SeqLike'
]


from . import strategies
from ._pydecoys import (
    from_fasta,
    to_fasta,
    from_SeqRecords,
    from_seqs,
    from_tuples,
    SeqRecord_as_decoy,
    seq_as_decoy,
    tuple_as_decoy,
    register,
    SeqLike
)

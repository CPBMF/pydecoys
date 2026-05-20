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
PyDecoys
========

Provides
    1. Functionality to generate decoys from proteins, including IO functions
    1. Integrated API to `Biopython`
    1. Easy implementation of custom decoy strategies within the API
    1. A CLI app to easily generate a decoy fasta file from a target fasta file

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
    Functions and API for decoy generation strategies

Avaliable decoy strategies
--------------------------
Each decoy strategy is specified by a lowercase string. The
:class:`strategies.DecoyGenerator` fns used for each strategy are available
under :mod:`strategies`. New strategies can be created by following the
:class:`strategies.DecoyGenerator` protocol or by instantiating
:class:`strategies.PseudoReverseRule` or :class:`strategies.PseudoShuffleRule`.

- reverse:                     Reverse protein
- reverse-keepn:               Reverse protein, except N-terminal aa
- reverse-keepc:               Reverse protein, except C-terminal aa
- reverse-keepterm:            Reserse protein, except terminal aas
- shuffle:                     Shuffle protein
- shuffle-keepn:               Shuffle protein, except N-terminal aa
- shuffle-keepc:               Shuffle protein, except C-terminal aa
- shuffle-keepterm:            Shuffle protein, except terminal aas
- pseudoreverse-trypsin:       Pseudo-reverse trypsin fragments
- pseudoreverse-stricttrypsin: Pseudo-reverse strict trypsin fragments
- pseudoreverse-argc:          Pseudo-reverse ArgC fragments
- pseudoreverse-aspn:          Pseudo-reverse AspN fragments
- pseudoreverse-chymo:         Pseudo-reverse chymotrypsin fragments
- pseudoreverse-gluc:          Pseudo-reverse GluC fragments
- pseudoreverse-lysc:          Pseudo-reverse LysC fragments
- pseudoreverse-lysn:          Pseudo-reverse LysN fragments
- pseudoreverse-stricttrypsin-keepn: Pseudo-reverse strict trypsin fragments, except N-terminal aa
- pseudoshuffle-trypsin:       Pseudo-shuffle trypsin fragments
- pseudoshuffle-stricttrypsin: Pseudo-shuffle strict trypsin fragments
- pseudoshuffle-argc:          Pseudo-shuffle ArgC fragments
- pseudoshuffle-aspn:          Pseudo-shuffle AspN fragments
- pseudoshuffle-chymo:         Pseudo-shuffle chymotrypsin fragments
- pseudoshuffle-gluc:          Pseudo-shuffle GluC fragments
- pseudoshuffle-lysc:          Pseudo-shuffle LysC fragments
- pseudoshuffle-lysn:          Pseudo-shuffle LysN fragments
- pseudoshuffle-stricttrypsin-keepn: Pseudo-shuffle strict trypsin fragments, except N-terminal aa
"""  # noqa: E501


__version_info__ = (0, 1, 0)
__version__ = '.'.join([str(i) for i in __version_info__])

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
    register
)

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

"""Module for handling decoy generation from target protein sequences.

This module exposes a series of functions to generate decoys from different
record representations, as well as related API. The IO functions require
`Biopython` to be installed.

- (fn) :func:`from_fasta`
- (fn) :func:`to_fasta`
- (fn) :func:`from_SeqRecords`
- (fn) :func:`from_seqs`
- (fn) :func:`from_tuples`
- (fn) :func:`SeqRecord_as_decoy`
- (fn) :func:`seq_as_decoy`
- (fn) :func:`tuple_as_decoy`
- (fn) :func:`register`
- (Protocol) :class:`DecoyGenerator`
- (TypeAlias) :obj:`SeqLike`
- (submodule) :mod:`DecoyStrategy`

Avaliable Decoy Strategies:
    Each decoy strategy is specified by a lowercase string. The
    :class:`DecoyGenerator` fns used for each strategy are available under
    :mod:`strategies`. New strategies can be created by following the
    :class:`DecoyGenerator` protocol or by instantiating
    :class:`PseudoReverseRule` or :class:`PseudoShuffleRule`.

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
    'SeqLike',
    'strategies',
]


from . import strategies
from .strategies import SeqLike
from ._decoys import (
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

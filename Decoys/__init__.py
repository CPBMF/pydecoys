# Copyright © 2026 Bruno Maestri A Becker
#
# This file is part of Decoys.
#
# Decoys is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Decoys is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# Decoys. If not, see <https://www.gnu.org/licenses/>.

"""Module for handling decoy generation from target protein sequences.

Attributes
----------
- (fn) :func:`generate`
- (fn) :func:`from_target`
- (fn) :func:`register`
- (TypeAlias) :obj:`DecoyGenerator`
- (submodule) :mod:`DecoyStrategy`

Avaliable Decoy Strategies
--------------------------
Each decoy strategy is specified by a lowercase string. The
:obj:`DecoyGenerator` fns used for each strategy are available under
:mod:`DecoyStrategy`.

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
- pseudoshuffle-trypsin:       Pseudo-shuffle trypsin fragments
- pseudoshuffle-stricttrypsin: Pseudo-shuffle strict trypsin fragments
- pseudoshuffle-argc:          Pseudo-shuffle ArgC fragments
- pseudoshuffle-aspn:          Pseudo-shuffle AspN fragments
- pseudoshuffle-chymo:         Pseudo-shuffle chymotrypsin fragments
- pseudoshuffle-gluc:          Pseudo-shuffle GluC fragments
- pseudoshuffle-lysc:          Pseudo-shuffle LysC fragments
- pseudoshuffle-lysn:          Pseudo-shuffle LysN fragments
"""


import typing as _t

from Bio.SeqRecord import SeqRecord

import DecoyStrategy
from .DecoyStrategy import DecoyGenerator
# Expose DecoyGenerator TypeAlias as part of the top-level API


__all__ = [
    'DecoyStrategy',
    'DecoyGenerator',
    'generate',
    'register',
]

__version_info__ = (0, 1, 0)
__version__ = '.'.join([str(i) for i in __version_info__])


_decoy_strategy: dict[str, DecoyGenerator] = {
    "reverse": DecoyStrategy.reverse,
    "reverse-keepn": DecoyStrategy.reverse_keep_n,
    "reverse-keepc": DecoyStrategy.reverse_keep_c,
    "reverse-keepterm": DecoyStrategy.reverse_keep_term,
    "shuffle": DecoyStrategy.shuffle,
    "shuffle-keepn": DecoyStrategy.shuffle_keep_n,
    "shuffle-keepc": DecoyStrategy.shuffle_keep_c,
    "shuffle-keepterm": DecoyStrategy.shuffle_keep_term,
    "pseudoreverse-trypsin": DecoyStrategy.pseudoreverse_trypsin,
    "pseudoreverse-stricttrypsin": DecoyStrategy.pseudoreverse_stricttrypsin,
    "pseudoreverse-argc": DecoyStrategy.pseudoreverse_argc,
    "pseudoreverse-aspn": DecoyStrategy.pseudoreverse_aspn,
    "pseudoreverse-chymo": DecoyStrategy.pseudoreverse_chymo,
    "pseudoreverse-gluc": DecoyStrategy.pseudoreverse_gluc,
    "pseudoreverse-lysc": DecoyStrategy.pseudoreverse_lysc,
    "pseudoreverse-lysn": DecoyStrategy.pseudoreverse_lysn,
    "pseudoshuffle-trypsin": DecoyStrategy.pseudoshuffle_trypsin,
    "pseudoshuffle-stricttrypsin": DecoyStrategy.pseudoshuffle_stricttrypsin,
    "pseudoshuffle-argc": DecoyStrategy.pseudoshuffle_argc,
    "pseudoshuffle-aspn": DecoyStrategy.pseudoshuffle_aspn,
    "pseudoshuffle-chymo": DecoyStrategy.pseudoshuffle_chymo,
    "pseudoshuffle-gluc": DecoyStrategy.pseudoshuffle_gluc,
    "pseudoshuffle-lysc": DecoyStrategy.pseudoshuffle_lysc,
    "pseudoshuffle-lysn": DecoyStrategy.pseudoshuffle_lysn,
}


def generate(
    sequences: _t.Iterable[SeqRecord] | SeqRecord,
    strategy: str,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> _t.Generator[SeqRecord, None, None]:
    """Lazily apply a decoy generation strategy to a set of sequences.

    Args:
        sequences: A list (or iterator) of :class:`Bio.SeqRecord.SeqRecord`
            objects, or a single :class:`Bio.SeqRecord.SeqRecord`.
        strategy: Lower case string especifying the decoy strategy to be used.
        decoy_tag: An optional tag that is to be appended to each input's
            :attr:`Bio.SeqRecord.SeqRecord.id`. Defaults to `'decoy_'`.
        prefix: If `False`, `decoy_tag` is suffixed, otherwise it's prefixed.
            Defaults to `True`.

    Yields:
        A decoy version of the next SeqRecord in `sequences`.

    Examples:
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Decoys import generate
        >>> seqs = [
        ...     SeqRecord('DNIDYKAVYR', 'seq1'),
        ...     SeqRecord('QSYMCTVTHP', 'seq2'),
        ...     SeqRecord('CQWSLTEELL', 'seq3'),
        ... ]
        >>> for decoy in generate(seqs, 'reverse'):
        ...     print(f'{decoy.id}: {decoy.seq}')
        decoy_seq1: Seq('RYVAKYDIND')
        decoy_seq2: Seq('PHTVTCMYSQ')
        decoy_seq3: Seq('LLEETLSWQC')
    """

    if not isinstance(strategy, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")

    if isinstance(sequences, SeqRecord):
        sequences = [sequences]

    decoy_generator = _decoy_strategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    for sequence in sequences:
        if sequence.seq is None:
            raise ValueError(f"Seq not present for SeqRecord '{sequence.id}'")

        id = sequence.id if sequence.id else ""
        id = decoy_tag + id if prefix else id + decoy_tag
        seq = decoy_generator(sequence.seq)

        yield SeqRecord(seq, id, description="")


def from_target(
    sequence: SeqRecord,
    strategy: str,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> SeqRecord:
    """Apply a decoy generation strategy to a given sequence.

    Args:
        sequence: A single :class:`Bio.SeqRecord.SeqRecord`.
        strategy: Lower case string especifying the decoy strategy to be used.
        decoy_tag: An optional tag that is to be appended to each input's
            :attr:`Bio.SeqRecord.SeqRecord.id`. Defaults to `'decoy_'`.
        prefix: If `False`, `decoy_tag` is suffixed, otherwise it's prefixed.
            Defaults to `True`.

    Returns:
        A decoy version of `sequence`.

    Examples:
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Decoys import from_target
        >>> seq = SeqRecord('DNIDYKAVYR', 'seq1')
        >>> decoy = from_target(seq, 'reverse')
        >>> print(f'{decoy.id}: {decoy.seq}')
        decoy_seq1: Seq('RYVAKYDIND')
    """

    if not isinstance(strategy, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")

    decoy_generator = _decoy_strategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    if sequence.seq is None:
        raise ValueError(f"Seq not present for SeqRecord '{sequence.id}'")

    id = sequence.id if sequence.id else ""
    id = decoy_tag + id if prefix else id + decoy_tag
    seq = decoy_generator(sequence.seq)

    return SeqRecord(seq, id, description="")


def register(strategy: str, decoy_generator_fn: DecoyGenerator) -> None:
    """Register a new decoy strategy that can be used with :func:`generate`.

    Args:
        strategy: Lower case string identifying the decoy strategy. Must not be
            already defined.
        decoy_generator_fn: A function that should take a :class:`Bio.Seq.Seq`
            object and return its decoy version.

    Examples:
        Given a :func:`random_seq` function that takes a Seq and returns a new,
        unrelated Seq:

        >>> from Bio.SeqRecord import SeqRecord
        >>> from Decoys import register, generate
        >>> register('randomseq', random_seq)
        >>> seq = SeqRecord('DNIDYKAVYR', 'seq1')
        >>> decoy = from_target(seq, 'randomseq')
        >>> print(f'{decoy.id}: {decoy.seq}')
        decoy_seq1: Seq('LLEETLSWQC')
    """

    if not isinstance(strategy, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    if strategy in _decoy_strategy:
        raise ValueError(f"Strategy '{strategy}' already exists")

    _decoy_strategy[strategy] = decoy_generator_fn

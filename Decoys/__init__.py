# Copyright © 2026 Bruno Maestri A Becker
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


import typing as _t

from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

from .DecoyStrategy import (
    reverse,
    reverse_keep_n,
    reverse_keep_c,
    reverse_keep_term,
    shuffle,
    PseudoReverseRule,
    PseudoShuffleRule,
    DecoyGenerator,
)


PseudoReverseTrypsin = PseudoReverseRule("KR", nocut="P")
PseudoReverseStrictTrypsin = PseudoReverseRule("KR")
PseudoReverseArgC = PseudoReverseRule("R", nocut="P")
PseudoReverseAspN = PseudoReverseRule("D", sense="N")
PseudoReverseChymo = PseudoReverseRule("FLWY", nocut="P")
PseudoReverseGluC = PseudoReverseRule("DE", nocut="P")
PseudoReverseLysC = PseudoReverseRule("K", nocut="P")
PseudoReverseLysN = PseudoReverseRule("K", sense="N")

PseudoShuffleTrypsin = PseudoShuffleRule("KR", nocut="P")
PseudoShuffleStrictTrypsin = PseudoShuffleRule("KR")
PseudoShuffleArgC = PseudoShuffleRule("R", nocut="P")
PseudoShuffleAspN = PseudoShuffleRule("D", sense="N")
PseudoShuffleChymo = PseudoShuffleRule("FLWY", nocut="P")
PseudoShuffleGluC = PseudoShuffleRule("DE", nocut="P")
PseudoShuffleLysC = PseudoShuffleRule("K", nocut="P")
PseudoShuffleLysN = PseudoShuffleRule("K", sense="N")


_DecoyStrategy: dict[str, DecoyGenerator] = {
    "reverse": reverse,
    "reverse-keepn": reverse_keep_n,
    "reverse-keepc": reverse_keep_c,
    "reverse-keepterm": reverse_keep_term,
    "shuffle": shuffle,
    "pseudoreverse-trypsin": PseudoReverseTrypsin,
    "pseudoreverse-stricttrypsin": PseudoReverseStrictTrypsin,
    "pseudoreverse-argc": PseudoReverseArgC,
    "pseudoreverse-aspn": PseudoReverseAspN,
    "pseudoreverse-chymo": PseudoReverseChymo,
    "pseudoreverse-gluc": PseudoReverseGluC,
    "pseudoreverse-lysc": PseudoReverseLysC,
    "pseudoreverse-lysn": PseudoReverseLysN,
    "pseudoshuffle-trypsin": PseudoShuffleTrypsin,
    "pseudoshuffle-stricttrypsin": PseudoShuffleStrictTrypsin,
    "pseudoshuffle-argc": PseudoReverseArgC,
    "pseudoshuffle-aspn": PseudoShuffleAspN,
    "pseudoshuffle-chymo": PseudoShuffleChymo,
    "pseudoshuffle-gluc": PseudoShuffleGluC,
    "pseudoshuffle-lysc": PseudoShuffleLysC,
    "pseudoshuffle-lysn": PseudoShuffleLysN,
}


def generate(
    sequences: _t.Iterable[SeqRecord] | SeqRecord,
    strategy: str,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> _t.Generator[SeqRecord, None, None]:
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

    decoy_generator = _DecoyStrategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    for sequence in sequences:
        if sequence.seq is None:
            raise ValueError(f"Seq not present for SeqRecord '{sequence.id}'")

        id = sequence.id if sequence.id else ""
        id = decoy_tag + id if prefix else id + decoy_tag
        seq = decoy_generator(sequence.seq)

        yield SeqRecord(seq, id, description="")


def register(strategy: str, fn: _t.Callable[[Seq | MutableSeq], Seq]) -> None:
    if not isinstance(strategy, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    if strategy in _DecoyStrategy:
        raise ValueError(f"Strategy '{strategy}' already exists")

    _DecoyStrategy[strategy] = fn

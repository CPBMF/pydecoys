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

"""Internal module responsible for lazily importing Biopython and registering
the relevant names.
"""

import re
from typing import Generator, Iterable

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq

import Decoys
from .DecoyStrategy import PseudoReverseRule, PseudoShuffleRule


def _decoys_from_SeqRecords(
    sequences: Iterable[SeqRecord] | SeqRecord,
    strategy: str,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> Generator[SeqRecord, None, None]:
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

    decoy_generator = Decoys._decoy_strategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    for sequence in sequences:
        if sequence.seq is None:
            raise ValueError(f"Seq not present for SeqRecord '{sequence.id}'")

        id = sequence.id if sequence.id else ""
        id = decoy_tag + id if prefix else id + decoy_tag
        seq = decoy_generator(sequence.seq)

        yield SeqRecord(seq, id, description="")


def _register():
    @PseudoReverseRule.__call__.register(Seq)         # type: ignore
    @PseudoReverseRule.__call__.register(MutableSeq)  # type: ignore
    def reverse_decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq:
        fragments = re.split(self._pattern, str(sequence))

        rev_frags = [frag[::-1] for frag in fragments]
        return Seq("".join(rev_frags))

    reverse_decoy_from_Seq.__doc__ = PseudoReverseRule.decoy_from_Seq.__doc__
    PseudoReverseRule.decoy_from_Seq = reverse_decoy_from_Seq

    @PseudoShuffleRule.__call__.register(Seq)         # type: ignore
    @PseudoShuffleRule.__call__.register(MutableSeq)  # type: ignore
    def shuffle_decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq:
        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return Seq("".join(shuf_frags))

    shuffle_decoy_from_Seq.__doc__ = PseudoShuffleRule.decoy_from_Seq.__doc__
    PseudoShuffleRule.decoy_from_Seq = shuffle_decoy_from_Seq

    _decoys_from_SeqRecords.__doc__ = Decoys.from_SeqRecords.__doc__
    Decoys.from_SeqRecords = _decoys_from_SeqRecords

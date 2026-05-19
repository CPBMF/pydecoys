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

"""Internal module responsible for lazily importing Biopython, registering
the relevant names and namespacing utility functions related to the Biopython
interface.
"""

import re
from typing import Iterable

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq

from .DecoyStrategy import PseudoReverseRule, PseudoShuffleRule


def SeqRecord_to_tuple(record: SeqRecord) -> tuple[str, str]:
    return (record.id if record.id else "", str(record.seq))


def tuple_to_SeqRecord(record: tuple[str, str]) -> SeqRecord:
    return SeqRecord(
        Seq(record[1]),
        record[0],
        name="",
        description=""
    )


def iter_SeqRecord(
    seq: Iterable[SeqRecord] | SeqRecord
) -> Iterable[SeqRecord]:
    if isinstance(seq, SeqRecord):
        return [seq]
    return seq


def str_to_Seq(sequence: str) -> Seq:
    return Seq(sequence)


def register():
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

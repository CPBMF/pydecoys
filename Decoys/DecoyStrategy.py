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


import random
import re
from typing import Callable, Literal, TypeAlias

from Bio.Seq import Seq, MutableSeq


DecoyGenerator: TypeAlias = Callable[[Seq | MutableSeq], Seq]


# So shuffled decoys are always reproducible
_rng = random.Random(10)


def reverse(sequence: Seq | MutableSeq) -> Seq:
    return sequence[::-1]


def reverse_keep_n(sequence: Seq | MutableSeq) -> Seq:
    return sequence[0] + sequence[:0:-1]


def reverse_keep_c(sequence: Seq | MutableSeq) -> Seq:
    return sequence[-2::-1] + sequence[-1]


def reverse_keep_term(sequence: Seq | MutableSeq) -> Seq:
    return sequence[0] + sequence[-2:0:-1] + sequence[-1]


def shuffle(sequence: Seq | MutableSeq) -> Seq:
    new = list(sequence)
    _rng.shuffle(new)
    return Seq("".join(new))


def shuffle_keep_n(sequence: Seq | MutableSeq) -> Seq:
    new = list(sequence[1:])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new))


def shuffle_keep_c(sequence: Seq | MutableSeq) -> Seq:
    new = list(sequence[:-1])
    _rng.shuffle(new)
    return Seq("".join(new) + sequence[-1])


def shuffle_keep_term(sequence: Seq | MutableSeq) -> Seq:
    new = list(sequence[1:-1])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new) + sequence[-1])


class PseudoReverseRule:
    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None
    ) -> None:
        if sense == 'N' and nocut is not None:
            raise ValueError("Cannot have nocut specification with sense N")

        pattern = rf"([{cut}])"
        if nocut is not None:
            pattern += rf"(?!{nocut})"
        self._pattern = re.compile(pattern)

    def __call__(self, sequence: Seq | MutableSeq) -> Seq:
        fragments = re.split(self._pattern, str(sequence))

        rev_frags = [frag[::-1] for frag in fragments]
        return Seq("".join(rev_frags))


class PseudoShuffleRule:
    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None
    ) -> None:
        if sense == 'N' and nocut is not None:
            raise ValueError("Cannot have nocut specification with sense N")

        pattern = rf"([{cut}])"
        if nocut is not None:
            pattern += rf"(?!{nocut})"
        self._pattern = re.compile(pattern)

    def __call__(self, sequence: Seq | MutableSeq) -> Seq:
        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return Seq("".join(shuf_frags))

    def _shuffle(self, frag: str) -> str:
        new = list(frag)
        _rng.shuffle(new)
        return "".join(new)

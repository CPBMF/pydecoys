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

"""Module that handles generating decoy proteins from target proteins.

This module serves as a namespace for all the :obj:`DecoyGenerator` functions.
"""


import random
import re
from collections.abc import Callable
from typing import Literal, TypeAlias

from Bio.Seq import Seq, MutableSeq


DecoyGenerator: TypeAlias = Callable[[Seq | MutableSeq], Seq]


# So shuffled decoys are always reproducible
_rng = random.Random(10)


def reverse(sequence: Seq | MutableSeq) -> Seq:
    """Return the reversed `sequence`."""
    return sequence[::-1]


def reverse_keep_n(sequence: Seq | MutableSeq) -> Seq:
    """Return the reversed `sequence`, except N-terminal."""
    return sequence[0] + sequence[:0:-1]


def reverse_keep_c(sequence: Seq | MutableSeq) -> Seq:
    """Return the reversed `sequence`, except C-terminal."""
    return sequence[-2::-1] + sequence[-1]


def reverse_keep_term(sequence: Seq | MutableSeq) -> Seq:
    """Return the reversed `sequence`, except terminals."""
    return sequence[0] + sequence[-2:0:-1] + sequence[-1]


def shuffle(sequence: Seq | MutableSeq) -> Seq:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    _rng.shuffle(new)
    return Seq("".join(new))


def shuffle_keep_n(sequence: Seq | MutableSeq) -> Seq:
    """Return the shuffled `sequence`, except N-terminal."""
    new = list(sequence[1:])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new))


def shuffle_keep_c(sequence: Seq | MutableSeq) -> Seq:
    """Return the shuffled `sequence`, except C-terminal."""
    new = list(sequence[:-1])
    _rng.shuffle(new)
    return Seq("".join(new) + sequence[-1])


def shuffle_keep_term(sequence: Seq | MutableSeq) -> Seq:
    """Return the shuffled `sequence`, except terminals."""
    new = list(sequence[1:-1])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new) + sequence[-1])


class PseudoReverseRule:
    """A :obj:`DecoyGenerator` that applies pseudo-reverse decoy generation
    with pre-specified enzyme specifications.

    A :class:`PseudoReverseRule` object is a callable object. The enzyme
    specifications can be checked with the attributes :attr:`cut`,
    :attr:`sense` and :attr:`nocut`.
    """

    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None
    ) -> None:
        """Initialize a pseudo-reverse :obj:`DecoyGenerator` with the specified
        enzyme specifications.

        Args:
            cut: Aminoacids the enzyme cleaves at as a string.
            sense: Sense of the enzyme (whether it cleaves the 'C' or 'N' bond
                of the cleavage site). Defaults to 'C'.
            nocut: Aminoacids that stop cleavage as a string, or `None`. If
                given, the enzyme won't cut aminoacids followed by these.
                Defaults to `None`.

        Examples:
            To specify a pseudo-reverse generator for trypsin:

            >>> from Bio.Seq import Seq
            >>> from Decoys import PseudoReverseRule
            >>> rev = PseudoReverseRule("KR", nocut="P")
            >>> print(rev.cut, rev.nocut, rev.sense, sep=', ')
            KR, P, C
        """

        if sense == 'N' and nocut is not None:
            raise ValueError("Cannot have nocut specification with sense N")

        self._cut = cut
        self._nocut = nocut
        # Without type hints this is cast to a str
        self._sense: Literal['C', 'N'] = sense

        pattern = rf"([{cut}])"
        if nocut is not None:
            pattern += rf"(?!{nocut})"
        self._pattern = re.compile(pattern)

    def __call__(self, sequence: Seq | MutableSeq) -> Seq:
        """Receive a :class:`Bio.Seq.Seq` and return a pseudo-reversed decoy.

        Args:
            sequence: A single :class:`Bio.Seq.Seq`.

        Returns:
            A pseudo-reversed version of `sequence`, according to the enzyme
                specifications given at class instantiation.

        Example:
            >>> from Bio.Seq import Seq
            >>> from Decoys import PseudoReverseRule
            >>> seq = Seq('QSYKPTRTHQ')
            >>> rev = PseudoReverseRule("KR", nocut="P")
            >>> print(rev(seq))
            Seq('TPKYSQRQHT')
        """

        fragments = re.split(self._pattern, str(sequence))

        rev_frags = [frag[::-1] for frag in fragments]
        return Seq("".join(rev_frags))

    @property
    def cut(self) -> str:
        """Cleavage sites as a string."""
        return self._cut

    @property
    def sense(self) -> Literal['C', 'N']:
        """Sense of cleavage."""
        return self._sense

    @property
    def nocut(self) -> str | None:
        """Aminoacids that stop cleavage as a string."""
        return self._nocut


class PseudoShuffleRule:
    """A :obj:`DecoyGenerator` that applies pseudo-shuffle decoy generation
    with pre-specified enzyme specifications.

    A :class:`PseudoShuffleRule` object is a callable object. The enzyme
    specifications can be checked with the attributes :attr:`cut`,
    :attr:`sense` and :attr:`nocut`.
    """

    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None
    ) -> None:
        """Initialize a pseudo-shuffle :obj:`DecoyGenerator` with the specified
        enzyme specifications.

        Args:
            cut: Aminoacids the enzyme cleaves at as a string.
            sense: Sense of the enzyme (whether it cleaves the 'C' or 'N' bond
                of the cleavage site). Defaults to 'C'.
            nocut: Aminoacids that stop cleavage as a string, or `None`. If
                given, the enzyme won't cut aminoacids followed by these.
                Defaults to `None`.

        Examples:
            To specify a pseudo-shuffle generator for trypsin:

            >>> from Bio.Seq import Seq
            >>> from Decoys import PseudoShuffleRule
            >>> shuf = PseudoShuffleRule("KR", nocut="P")
            >>> print(shuf.cut, shuf.nocut, shuf.sense, sep=', ')
            KR, P, C
        """

        if sense == 'N' and nocut is not None:
            raise ValueError("Cannot have nocut specification with sense N")

        self._cut = cut
        self._nocut = nocut
        # Without type hint sense is cast to str
        self._sense: Literal['C', 'N'] = sense

        pattern = rf"([{cut}])"
        if nocut is not None:
            pattern += rf"(?!{nocut})"
        self._pattern = re.compile(pattern)

    def __call__(self, sequence: Seq | MutableSeq) -> Seq:
        """Receive a :class:`Bio.Seq.Seq` and return a pseudo-shuffled decoy.

        Args:
            sequence: A single :class:`Bio.Seq.Seq`.

        Returns:
            A pseudo-reversed version of `sequence`, according to the enzyme
                specifications given at class instantiation.

        Example:
            >>> from Bio.Seq import Seq
            >>> from Decoys import PseudoShuffleRule
            >>> shuf = Seq('QSYKPTRTHQ')
            >>> shuf = PseudoShuffleRule("KR", nocut="P")
            >>> print(shuf(seq))
            Seq('YTSKQPRQHT')
        """

        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return Seq("".join(shuf_frags))

    @property
    def cut(self) -> str:
        """Cleavage sites as a string."""
        return self._cut

    @property
    def sense(self) -> Literal['C', 'N']:
        """Sense of cleavage."""
        return self._sense

    @property
    def nocut(self) -> str | None:
        """Aminoacids that stop cleavage as a string."""
        return self._nocut

    def _shuffle(self, frag: str) -> str:
        new = list(frag)
        _rng.shuffle(new)
        return "".join(new)


# Pre-defined pseudo-reverse and pseudo-shuffle DecoyGenerators
pseudoreverse_trypsin: DecoyGenerator = PseudoReverseRule("KR", nocut="P")
pseudoreverse_stricttrypsin: DecoyGenerator = PseudoReverseRule("KR")
pseudoreverse_argc: DecoyGenerator = PseudoReverseRule("R", nocut="P")
pseudoreverse_aspn: DecoyGenerator = PseudoReverseRule("D", sense="N")
pseudoreverse_chymo: DecoyGenerator = PseudoReverseRule("FLWY", nocut="P")
pseudoreverse_gluc: DecoyGenerator = PseudoReverseRule("DE", nocut="P")
pseudoreverse_lysc: DecoyGenerator = PseudoReverseRule("K", nocut="P")
pseudoreverse_lysn: DecoyGenerator = PseudoReverseRule("K", sense="N")

pseudoshuffle_trypsin: DecoyGenerator = PseudoShuffleRule("KR", nocut="P")
pseudoshuffle_stricttrypsin: DecoyGenerator = PseudoShuffleRule("KR")
pseudoshuffle_argc: DecoyGenerator = PseudoShuffleRule("R", nocut="P")
pseudoshuffle_aspn: DecoyGenerator = PseudoShuffleRule("D", sense="N")
pseudoshuffle_chymo: DecoyGenerator = PseudoShuffleRule("FLWY", nocut="P")
pseudoshuffle_gluc: DecoyGenerator = PseudoShuffleRule("DE", nocut="P")
pseudoshuffle_lysc: DecoyGenerator = PseudoShuffleRule("K", nocut="P")
pseudoshuffle_lysn: DecoyGenerator = PseudoShuffleRule("K", sense="N")

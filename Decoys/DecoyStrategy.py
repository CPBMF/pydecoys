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

from __future__ import annotations

from functools import singledispatchmethod
import random
import re
from typing import Literal, Protocol, TypeAlias, TYPE_CHECKING, overload

if TYPE_CHECKING:
    from Bio.Seq import Seq, MutableSeq


SeqLike: TypeAlias = 'str | Seq | MutableSeq'


class DecoyGenerator(Protocol):
    """Protocol defining a decoy generator function that applies a decoy
    strategy.

    Classes following this protocol should implement `__call__` as a
    `Callable[[SeqLike], SeqLile]`. For most cases, a :obj:`SeqLike` object
    can be treated as a `str`, including splicing. Regardless, for concistency
    it is recommended that `__call__` be overloaded as follows:

    - `Callable[[str], str]`
    - `Callable[[Seq | MutableSeq], Seq]` or
    `Callable[[Seq | MutableSeq], Seq | MutableSeq]`
    """

    @overload
    def __call__(self, sequence: Seq | MutableSeq) -> Seq:
        ...

    @overload
    def __call__(self, sequence: str) -> str:
        ...

    def __call__(self, sequence: SeqLike) -> SeqLike:
        """Generate a decoy version of a given sequence.

        Args:
            sequence: A target sequence.

        Returns:
            The decoy version of `sequence`.
        """
        ...


# So shuffled decoys are always reproducible
_rng = random.Random(10)


def reverse(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`."""
    return sequence[::-1]


def reverse_keep_n(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except N-terminal."""
    return sequence[0] + sequence[:0:-1]


def reverse_keep_c(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except C-terminal."""
    return sequence[-2::-1] + sequence[-1]


def reverse_keep_term(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except terminals."""
    return sequence[0] + sequence[-2:0:-1] + sequence[-1]


def shuffle(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    _rng.shuffle(new)
    return Seq("".join(new))


def shuffle_keep_n(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except N-terminal."""
    new = list(sequence[1:])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new))


def shuffle_keep_c(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except C-terminal."""
    new = list(sequence[:-1])
    _rng.shuffle(new)
    return Seq("".join(new) + sequence[-1])


def shuffle_keep_term(sequence: SeqLike) -> SeqLike:
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

    @singledispatchmethod
    def __call__(self, sequence: SeqLike) -> SeqLike:
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
        from . import _Bio
        _Bio._register()
        return self.__call__(sequence)

    @__call__.register
    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to `PseudoReverseRule(sequence)`
        where `sequence` is a `str`, but avoids `singledispatch` overhead.
        """
        fragments = re.split(self._pattern, str(sequence))

        rev_frags = [frag[::-1] for frag in fragments]
        return "".join(rev_frags)

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq:
        """Convenience funcion. Equivalent to `PseudoReverseRule(sequence)`
        where `sequence` is a :class:`Bio.Seq.Seq`, but avoids `singledispatch`
        overhead.
        """
        from . import _Bio
        _Bio._register()
        return self.decoy_from_Seq(sequence)

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

    @singledispatchmethod
    def __call__(self, sequence: SeqLike) -> SeqLike:
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
        from . import _Bio
        _Bio._register()
        return self.__call__(sequence)

    @__call__.register
    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to `PseudoShuffleRule(sequence)`
        where `sequence` is a `str`, but avoids `singledispatch` overhead.
        """
        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return "".join(shuf_frags)

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq:
        """Convenience funcion. Equivalent to `PseudoShuffleRule(sequence)`
        where `sequence` is a :class:`Bio.Seq.Seq`, but avoids `singledispatch`
        overhead.
        """
        from . import _Bio
        _Bio._register()
        return self.decoy_from_Seq(sequence)

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
pseudoreverse_trypsin = PseudoReverseRule("KR", nocut="P")  # type: ignore # noqa: E501
pseudoreverse_stricttrypsin = PseudoReverseRule("KR")       # type: ignore # noqa: E501
pseudoreverse_argc = PseudoReverseRule("R", nocut="P")      # type: ignore # noqa: E501
pseudoreverse_aspn = PseudoReverseRule("D", sense="N")      # type: ignore # noqa: E501
pseudoreverse_chymo = PseudoReverseRule("FLWY", nocut="P")  # type: ignore # noqa: E501
pseudoreverse_gluc = PseudoReverseRule("DE", nocut="P")     # type: ignore # noqa: E501
pseudoreverse_lysc = PseudoReverseRule("K", nocut="P")      # type: ignore # noqa: E501
pseudoreverse_lysn = PseudoReverseRule("K", sense="N")      # type: ignore # noqa: E501

pseudoshuffle_trypsin = PseudoShuffleRule("KR", nocut="P")  # type: ignore # noqa: E501
pseudoshuffle_stricttrypsin = PseudoShuffleRule("KR")       # type: ignore # noqa: E501
pseudoshuffle_argc = PseudoShuffleRule("R", nocut="P")      # type: ignore # noqa: E501
pseudoshuffle_aspn = PseudoShuffleRule("D", sense="N")      # type: ignore # noqa: E501
pseudoshuffle_chymo = PseudoShuffleRule("FLWY", nocut="P")  # type: ignore # noqa: E501
pseudoshuffle_gluc = PseudoShuffleRule("DE", nocut="P")     # type: ignore # noqa: E501
pseudoshuffle_lysc = PseudoShuffleRule("K", nocut="P")      # type: ignore # noqa: E501
pseudoshuffle_lysn = PseudoShuffleRule("K", sense="N")      # type: ignore # noqa: E501

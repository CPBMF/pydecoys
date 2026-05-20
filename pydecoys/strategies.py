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
strategies
==========

Decoy-generation logic, as well as API to write new decoy strategies.

The main API of `strategies` is the :class:`DecoyGenerator` type. This is a
simple protocol that only implements a `__call__` function and appropriate
type overloads. For decoy strategies that need context from the target
database (for example, that use a Markov State Model), implement the
:class:`ContextfulGenerator` protocol. The :class:`PseudoReverseRule` and
:class:`PseudoShuffleRule` classes allow easy definition of new enzyme
specifications for pseudo-reverse and pseudo-shuffle strategies via
instantiation.

Available enzymes
-----------------
A list of pre-instantiated pseudo-reverse and pseudo-shuffle generators
that are available, following the name scheme `pseudoreverse_<enzyme>` and
`pseudoshuffle_<enzyme>`.

- pseudoreverse_trypsin:             cut KR, nocut P, sense C, keep_n False
- pseudoreverse_stricttrypsin:       cut KR, nocut None, sense C, keep_n False
- pseudoreverse_argc:                cut R, nocut P, sense C, keep_n False
- pseudoreverse_aspn:                cut D, nocut None, sense N, keep_n False
- pseudoreverse_chymo:               cut FLWY, nocut P, sense C, keep_n False
- pseudoreverse_gluc:                cut DE, nocut P, sense C, keep_n False
- pseudoreverse_lysc:                cut K, nocut P, sense C, keep_n False
- pseudoreverse_lysn:                cut K, nocut None, sense N
- pseudoreverse_stricttrypsin_keepn: cut KR, nocut None, sense C, keepn True
- pseudoshuffle_trypsin:             cut KR, nocut P, sense C, keep_n False
- pseudoshuffle_stricttrypsin:       cut KR, nocut None, sense C, keep_n False
- pseudoshuffle_argc:                cut R, nocut P, sense C, keep_n False
- pseudoshuffle_aspn:                cut D, nocut None, sense N, keep_n False
- pseudoshuffle_chymo:               cut FLWY, nocut P, sense C, keep_n False
- pseudoshuffle_gluc:                cut DE, nocut P, sense C, keep_n False
- pseudoshuffle_lysc:                cut K, nocut P, sense C, keep_n False
- pseudoshuffle_lysn:                cut K, nocut None, sense N
- pseudoshuffle_stricttrypsin_keepn: cut KR, nocut None, sense C, keepn True
"""

from __future__ import annotations

from functools import singledispatchmethod
import random
import re
from typing import (
    Literal,
    Protocol,
    TYPE_CHECKING,
    overload,
    runtime_checkable
)

if TYPE_CHECKING:
    from Bio.Seq import Seq, MutableSeq


type SeqLike = 'str | Seq | MutableSeq'
"""`SeqLike` objects can be indexed and spliced; `str` at runtime."""


class DecoyGenerator(Protocol):
    """Protocol defining a decoy generator function that applies a decoy
    strategy.

    Classes following this protocol should implement `__call__` as a
    `Callable[[SeqLike], SeqLike]`. For most cases, a :obj:`SeqLike` object
    can be treated as a `str`, including splicing and concatenation:
    ``seq[1::-1]`` and ``seq1 + seq2``.
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


@runtime_checkable
class ContextfulGenerator(DecoyGenerator, Protocol):
    """Protocol defining a decoy generator function that uses previously
    learned context.
    """

    def learn_context(
        self,
        sequences: list[SeqLike]
    ) -> None:
        """Receive the target proteins set to generate the necessary context.

        Parameters
        ----------
        sequences
            The target dataset.
        """
        ...


# So shuffled decoys are always reproducible
_rng = random.Random(10)


def reverse(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`."""
    return sequence[::-1]


def reverse_keep_n(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except N-terminal aa."""
    return sequence[0] + sequence[:0:-1]


def reverse_keep_c(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except C-terminal aa."""
    return sequence[-2::-1] + sequence[-1]


def reverse_keep_term(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except terminal aas."""
    return sequence[0] + sequence[-2:0:-1] + sequence[-1]


def shuffle(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    _rng.shuffle(new)
    return Seq("".join(new))


def shuffle_keep_n(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except N-terminal aa."""
    new = list(sequence[1:])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new))


def shuffle_keep_c(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except C-terminal aa."""
    new = list(sequence[:-1])
    _rng.shuffle(new)
    return Seq("".join(new) + sequence[-1])


def shuffle_keep_term(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except terminal aas."""
    new = list(sequence[1:-1])
    _rng.shuffle(new)
    return Seq(sequence[0] + "".join(new) + sequence[-1])


class PseudoReverseRule:
    """Appliy pseudo-reverse decoy generation with the specified enzyme
    properties.

    Callable object. Enzyme specifications can be checked via its attributes.

    Parameters
    ----------
    cut
        Cleavage sites as a string.
    sense
        Sense cleavage (whether the enzyme cleaves the 'C' or 'N' bond of the
        cleavage site).
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids followed by these..
    keep_n
        If `True`, the N-terminal aa isn't reverted.

    Examples
    --------
    >>> from pydecoys.strategies import PseudoReverseRule
    >>> rev = PseudoReverseRule("KR", nocut="P")
    >>> print(rev.cut, rev.nocut, rev.sense, sep=', ')
    KR, P, C
    >>> rev('QSYKPTRTHQ')
    'TPKYSQRQHT'
    """

    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None,
        keep_n: bool = False
    ) -> None:
        if sense == 'N' and nocut is not None:
            raise ValueError("Cannot have nocut specification with N sense")

        self._cut = cut
        self._nocut = nocut
        # Without type hints this is cast to a str
        self._sense: Literal['C', 'N'] = sense
        self._keep_n: bool = keep_n

        pattern = rf"([{cut}])"
        if self._keep_n:
            pattern = f"(^.|[{cut}])"
        if nocut is not None:
            pattern += rf"(?!{nocut})"
        self._pattern = re.compile(pattern)

    @singledispatchmethod
    def __call__(self, sequence: SeqLike) -> SeqLike:
        """Receive a sequence and return a pseudo-reversed decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        A pseudo-reversed version of `sequence`, according to the enzyme
        specifications given at class instantiation.

        Examples
        --------
        >>> from pydecoys.strategies import PseudoReverseRule
        >>> rev = PseudoReverseRule("KR", nocut="P")
        >>> print(rev.cut, rev.nocut, rev.sense, sep=', ')
        KR, P, C
        >>> rev('QSYKPTRTHQ')
        'TPKYSQRQHT'
        """
        from . import _bio
        _bio.register()
        return self.__call__(sequence)

    @__call__.register
    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to `PseudoReverseRule(sequence)`
        where `sequence` is a `str`, but avoids
        :class:`typing.singledispatchmethod` overhead.
        """
        fragments = re.split(self._pattern, str(sequence))

        rev_frags = [frag[::-1] for frag in fragments]
        return "".join(rev_frags)

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq:
        """Convenience funcion. Equivalent to `PseudoReverseRule(sequence)`
        where `sequence` is a `Seq`,  but avoids
        :class:`typing.singledispatchmethod` overhead.
        """
        from . import _bio
        _bio.register()
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

    @property
    def keep_n(self) -> bool:
        """If `True`, the N-terminal aa isn't reverted."""
        return self._keep_n


class PseudoShuffleRule:
    """Appliy pseudo-shuffle decoy generation with the specified enzyme
    properties.

    Callable object. Enzyme specifications can be checked via its attributes.

    Parameters
    ----------
    cut
        Cleavage sites as a string.
    sense
        Sense cleavage (whether the enzyme cleaves the 'C' or 'N' bond of the
        cleavage site).
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids followed by these..
    keep_n
        If `True`, the N-terminal aa isn't reverted.

    Examples
    --------
    >>> from pydecoys.strategies import PseudoShuffleRule
    >>> shuf = PseudoShuffleRule("KR", nocut="P")
    >>> print(shuf.cut, shuf.nocut, shuf.sense, sep=', ')
    KR, P, C
    >>> shuf('QSYKPTRTHQ')
    'YTSKQPRQHT'
    """

    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None,
        keep_n: bool = False,
    ) -> None:
        if sense == 'N' and nocut is not None:
            raise ValueError("Cannot have nocut specification with sense N")

        self._cut = cut
        self._nocut = nocut
        # Without type hints this is cast to a str
        self._sense: Literal['C', 'N'] = sense
        self._keep_n: bool = keep_n

        pattern = rf"([{cut}])"
        if self._keep_n:
            pattern = f"(^.|[{cut}])"
        if nocut is not None:
            pattern += rf"(?!{nocut})"
        self._pattern = re.compile(pattern)

    @singledispatchmethod
    def __call__(self, sequence: SeqLike) -> SeqLike:
        """Receive a sequence and return a pseudo-shuffled decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        A pseudo-shuffled version of `sequence`, according to the enzyme
        specifications given at class instantiation.

        Examples
        --------
        >>> from pydecoys.strategies import PseudoShuffleRule
        >>> shuf = PseudoShuffleRule("KR", nocut="P")
        >>> print(shuf.cut, shuf.nocut, shuf.sense, sep=', ')
        KR, P, C
        >>> shuf('QSYKPTRTHQ')
        'YTSKQPRQHT'
        """
        from . import _bio
        _bio.register()
        return self.__call__(sequence)

    @__call__.register
    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to `PseudoShuffleRule(sequence)`
        where `sequence` is a `str`,  but avoids
        :class:`typing.singledispatchmethod` overhead.
        """
        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return "".join(shuf_frags)

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq:
        """Convenience funcion. Equivalent to `PseudoShuffleRule(sequence)`
        where `sequence` is a :class:`Bio.Seq.Seq`,  but avoids
        :class:`typing.singledispatchmethod` overhead.
        """
        from . import _bio
        _bio.register()
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

    @property
    def keep_n(self) -> bool:
        """If `True`, the N-terminal aa isn't reverted."""
        return self._keep_n

    def _shuffle(self, frag: str) -> str:
        new = list(frag)
        _rng.shuffle(new)
        return "".join(new)


# Pre-defined pseudo-reverse and pseudo-shuffle DecoyGenerators
pseudoreverse_trypsin = PseudoReverseRule("KR", nocut="P")
pseudoreverse_stricttrypsin = PseudoReverseRule("KR")
pseudoreverse_argc = PseudoReverseRule("R", nocut="P")
pseudoreverse_aspn = PseudoReverseRule("D", sense="N")
pseudoreverse_chymo = PseudoReverseRule("FLWY", nocut="P")
pseudoreverse_gluc = PseudoReverseRule("DE", nocut="P")
pseudoreverse_lysc = PseudoReverseRule("K", nocut="P")
pseudoreverse_lysn = PseudoReverseRule("K", sense="N")
pseudoreverse_stricttrypsin_keepn = PseudoReverseRule("KR", keep_n=True)

pseudoshuffle_trypsin = PseudoShuffleRule("KR", nocut="P")
pseudoshuffle_stricttrypsin = PseudoShuffleRule("KR")
pseudoshuffle_argc = PseudoShuffleRule("R", nocut="P")
pseudoshuffle_aspn = PseudoShuffleRule("D", sense="N")
pseudoshuffle_chymo = PseudoShuffleRule("FLWY", nocut="P")
pseudoshuffle_gluc = PseudoShuffleRule("DE", nocut="P")
pseudoshuffle_lysc = PseudoShuffleRule("K", nocut="P")
pseudoshuffle_lysn = PseudoShuffleRule("K", sense="N")
pseudoshuffle_stricttrypsin_keepn = PseudoShuffleRule("KR", keep_n=True)

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

"""Decoy-generation logic, as well as API to write new decoy strategies.

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
Pre-initialized pseudo-reversers and pseudo-shufflers covering most proteases
are available, following the name scheme `pseudoreverse_<enzyme>` and
`pseudoshuffle_<enzyme>`.
"""

from __future__ import annotations

from collections.abc import Sequence
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


# This file has a lot of overloads, but sphinx autodoc doesn't seem to catch
# the overloads if I move them to a .pyi file.


type SeqLike = 'str | Seq | MutableSeq'
"""`SeqLike` objects can be indexed and spliced; `str` at runtime."""

type Seq_ = 'Seq'
"""`Seq` type that doesn't require Biopython; `str` at runtime."""

type MutableSeq_ = 'MutableSeq'
"""`MutableSeq` type that doesn't require Biopython; `str` at runtime."""


class DecoyGenerator(Protocol):
    """Protocol defining a decoy generator function that applies a decoy
    strategy.

    Classes following this protocol should implement `__call__` as a
    `Callable[[SeqLike], SeqLike]`. For most cases, a :obj:`SeqLike` object
    can be treated as a `str`, including splicing and concatenation:
    ``seq[1::-1]`` and ``seq1 + seq2``.
    """

    @overload
    def __call__(self, sequence: Seq) -> Seq:
        ...

    @overload
    def __call__(self, sequence:  MutableSeq) -> MutableSeq:
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
        sequences: Sequence[SeqLike]
    ) -> None:
        """Receive the target proteins set to generate the necessary context.

        Parameters
        ----------
        sequences
            The target dataset.
        """
        ...


# So shuffled decoys are always reproducible
rand = random.Random(10)
"""Random number generator for stochastic decoy strategies."""


@overload
def reverse(sequence: Seq) -> Seq: ...


@overload
def reverse(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse(sequence: str) -> str: ...


def reverse(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`."""
    return sequence[::-1]


@overload
def reverse_keep_n(sequence: Seq) -> Seq: ...


@overload
def reverse_keep_n(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse_keep_n(sequence: str) -> str: ...


def reverse_keep_n(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except N-terminal aa."""
    return sequence[0] + sequence[:0:-1]


@overload
def reverse_keep_c(sequence: Seq) -> Seq: ...


@overload
def reverse_keep_c(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse_keep_c(sequence: str) -> str: ...


def reverse_keep_c(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except C-terminal aa."""
    return sequence[-2::-1] + sequence[-1]


@overload
def reverse_keep_term(sequence: Seq) -> Seq: ...


@overload
def reverse_keep_term(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse_keep_term(sequence: str) -> str: ...


def reverse_keep_term(sequence: SeqLike) -> SeqLike:
    """Return the reversed `sequence`, except terminal aas."""
    return sequence[0] + sequence[-2:0:-1] + sequence[-1]


@overload
def shuffle(sequence: Seq) -> Seq: ...


@overload
def shuffle(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle(sequence: str) -> str: ...


def shuffle(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    rand.shuffle(new)
    return _cls_cast(sequence, "".join(new))


@overload
def shuffle_keep_n(sequence: Seq) -> Seq: ...


@overload
def shuffle_keep_n(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle_keep_n(sequence: str) -> str: ...


def shuffle_keep_n(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except N-terminal aa."""
    new = list(sequence[1:])
    rand.shuffle(new)
    return _cls_cast(sequence, sequence[0] + "".join(new))


@overload
def shuffle_keep_c(sequence: Seq) -> Seq: ...


@overload
def shuffle_keep_c(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle_keep_c(sequence: str) -> str: ...


def shuffle_keep_c(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except C-terminal aa."""
    new = list(sequence[:-1])
    rand.shuffle(new)
    return _cls_cast(sequence, "".join(new) + sequence[-1])


@overload
def shuffle_keep_term(sequence: Seq) -> Seq: ...


@overload
def shuffle_keep_term(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle_keep_term(sequence: str) -> str: ...


def shuffle_keep_term(sequence: SeqLike) -> SeqLike:
    """Return the shuffled `sequence`, except terminal aas."""
    new = list(sequence[1:-1])
    rand.shuffle(new)
    return _cls_cast(sequence, sequence[0] + "".join(new) + sequence[-1])


# Hackish solution, but it allows the code to always return the correct type
# without importing Biopython or deferring to another module
def _cls_cast(obj: SeqLike, sequence: str) -> SeqLike:
    cls = type(obj)
    return cls(sequence)


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
        enzyme won't cut aminoacids followed by these.
    keep_n
        If `True`, the N-terminal aa isn't reverted.

    Examples
    --------
    >>> from pydecoys.strategies import PseudoReverseRule
    >>> rev = PseudoReverseRule("KR", nocut="P")
    >>> print(rev.cut, rev.nocut, rev.sense, sep=', ')
    KR, P, C
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

    @overload
    def __call__(self, sequence: Seq) -> Seq: ...

    @overload
    def __call__(self, sequence: MutableSeq) -> MutableSeq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

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
        >>> rev('QSYKPTRTHQ')
        'TPKYSQRQHT'
        """
        fragments = re.split(self._pattern, str(sequence))
        rev_frags = [frag[::-1] for frag in fragments]
        return _cls_cast(sequence, "".join(rev_frags))

    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to ``PseudoReverseRule(sequence)``
        where `sequence` is a `str`.
        """
        fragments = re.split(self._pattern, str(sequence))
        rev_frags = [frag[::-1] for frag in fragments]
        return "".join(rev_frags)

    @overload
    def decoy_from_Seq(self, sequence: Seq) -> Seq: ...

    @overload
    def decoy_from_Seq(self, sequence: MutableSeq) -> MutableSeq: ...

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq | MutableSeq:
        """Convenience funcion. Equivalent to ``PseudoReverseRule(sequence)``
        where `sequence` is a `Seq` or `MutableSeq`.
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
        enzyme won't cut aminoacids followed by these.
    keep_n
        If `True`, the N-terminal aa isn't reverted.

    Examples
    --------
    >>> from pydecoys.strategies import PseudoShuffleRule
    >>> shuf = PseudoShuffleRule("KR", nocut="P")
    >>> print(shuf.cut, shuf.nocut, shuf.sense, sep=', ')
    KR, P, C
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

    @overload
    def __call__(self, sequence: Seq) -> Seq: ...

    @overload
    def __call__(self, sequence: MutableSeq) -> MutableSeq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

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
        >>> shuf('QSYKPTRTHQ')
        'YTSKQPRQHT'
        """
        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return _cls_cast(sequence, "".join(shuf_frags))

    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to ``PseudoShuffleRule(sequence)``
        where `sequence` is a `str`.
        """
        fragments = re.split(self._pattern, str(sequence))
        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return "".join(shuf_frags)

    @overload
    def decoy_from_Seq(self, sequence: Seq) -> Seq: ...

    @overload
    def decoy_from_Seq(self, sequence: MutableSeq) -> MutableSeq: ...

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq | MutableSeq:
        """Convenience funcion. Equivalent to ``PseudoShuffleRule(sequence)``
        where `sequence` is a `Seq` or `MutableSeq`.
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
        rand.shuffle(new)
        return "".join(new)


# Pre-defined pseudo-reverse and pseudo-shuffle DecoyGenerators
pseudoreverse_trypsin: DecoyGenerator = PseudoReverseRule("KR", nocut="P")
"""Enzymatic specifications - ``cut='KR'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""    # noqa: E501

pseudoreverse_stricttrypsin: DecoyGenerator = PseudoReverseRule("KR")
"""Enzymatic specifications - ``cut='KR'``, ``nocut=None``, ``sense='C'``; ``keep_n=False``."""   # noqa: E501

pseudoreverse_argc: DecoyGenerator = PseudoReverseRule("R", nocut="P")
"""Enzymatic specifications - ``cut='R'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""     # noqa: E501

pseudoreverse_aspn: DecoyGenerator = PseudoReverseRule("D", sense="N")
"""Enzymatic specifications - ``cut='D'``, ``nocut=None``, ``sense='N'``; ``keep_n=False``."""    # noqa: E501

pseudoreverse_chymo: DecoyGenerator = PseudoReverseRule("FLWY", nocut="P")
"""Enzymatic specifications - ``cut='FLWY'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""  # noqa: E501

pseudoreverse_gluc: DecoyGenerator = PseudoReverseRule("DE", nocut="P")
"""Enzymatic specifications - ``cut='DE'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""    # noqa: E501

pseudoreverse_lysc: DecoyGenerator = PseudoReverseRule("K", nocut="P")
"""Enzymatic specifications - ``cut='K'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""     # noqa: E501

pseudoreverse_lysn: DecoyGenerator = PseudoReverseRule("K", sense="N")
"""Enzymatic specifications - ``cut='K'``, ``nocut=None``, ``sense='N'``; ``keep_n=False``."""    # noqa: E501

pseudoreverse_stricttrypsin_keepn: DecoyGenerator = PseudoReverseRule("KR", keep_n=True)          # noqa: E501
"""Enzymatic specifications - ``cut='KR'``, ``nocut=None``, ``sense='C'``; ``keep_n=True``."""    # noqa: E501

pseudoshuffle_trypsin: DecoyGenerator = PseudoShuffleRule("KR", nocut="P")
"""Enzymatic specifications - ``cut='KR'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""    # noqa: E501

pseudoshuffle_stricttrypsin: DecoyGenerator = PseudoShuffleRule("KR")
"""Enzymatic specifications - ``cut='KR'``, ``nocut=None``, ``sense='C'``; ``keep_n=False``."""   # noqa: E501

pseudoshuffle_argc: DecoyGenerator = PseudoShuffleRule("R", nocut="P")
"""Enzymatic specifications - ``cut='R'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""     # noqa: E501

pseudoshuffle_aspn: DecoyGenerator = PseudoShuffleRule("D", sense="N")
"""Enzymatic specifications - ``cut='D'``, ``nocut=None``, ``sense='N'``; ``keep_n=False``."""    # noqa: E501

pseudoshuffle_chymo: DecoyGenerator = PseudoShuffleRule("FLWY", nocut="P")
"""Enzymatic specifications - ``cut='FLWY'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""  # noqa: E501

pseudoshuffle_gluc: DecoyGenerator = PseudoShuffleRule("DE", nocut="P")
"""Enzymatic specifications - ``cut='DE'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""    # noqa: E501

pseudoshuffle_lysc: DecoyGenerator = PseudoShuffleRule("K", nocut="P")
"""Enzymatic specifications - ``cut='K'``, ``nocut='P'``, ``sense='C'``; ``keep_n=False``."""     # noqa: E501

pseudoshuffle_lysn: DecoyGenerator = PseudoShuffleRule("K", sense="N")
"""Enzymatic specifications - ``cut='K'``, ``nocut=None``, ``sense='N'``; ``keep_n=False``."""    # noqa: E501

pseudoshuffle_stricttrypsin_keepn: DecoyGenerator = PseudoShuffleRule("KR", keep_n=True)          # noqa: E501
"""Enzymatic specifications - ``cut='KR'``, ``nocut=None``, ``sense='C'``; ``keep_n=True``."""    # noqa: E501

# Copyright (C) 2026 CPBMF and INCT-TB, PUCRS, Porto Alegre, Brazil
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

"""Core API (protocols, ABCs, types and utility) for :mod:`strategies`."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
import random
import re
from typing import (
    Final,
    Literal,
    Protocol,
    overload,
    TYPE_CHECKING,
    override,
    runtime_checkable,
)

if TYPE_CHECKING:
    from Bio.Seq import Seq, MutableSeq


type SeqLike = 'str | Seq | MutableSeq'
"""`SeqLike` objects can be indexed and spliced; `str` at runtime."""

type Seq_ = 'Seq'
"""`Seq` type that doesn't require Biopython; `str` at runtime."""

type MutableSeq_ = 'MutableSeq'
"""`MutableSeq` type that doesn't require Biopython; `str` at runtime."""


# So shuffled decoys are always reproducible
RAND: Final = random.Random(10)
"""Random number generator for stochastic decoy strategies."""

AMINOACIDS: Final = 'QWERTYIPASDFGHKLCVNM'
"""Standard 20 aminoacids single-letter codes, majuscule."""


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

    def learn_context(self, sequences: Sequence[SeqLike]) -> None:
        """Receive the target proteins set to generate the necessary context.

        Parameters
        ----------
        sequences
            The target dataset.
        """
        ...


class EnzymeSpecificGenerator(ABC):
    """Abstract base class for enzymatic aware decoy generation.

    This class creates a compiled regex pattern at instantiation that captures
    peptides that shouldn't be altered (cleavage sites, and maybe N- and
    C-termini if specified at class instantiation). The pattern can be accessed
    through ``self._pattern``.

    This class also save the enzymatic specifications as get-only attributes.

    Parameters
    ----------
    cut
        Cleavage sites as a string.
    sense
        Whether the enzyme cleaves the C or N bond of the cleavage site.
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids with a C-terminal followed by these.
    keep_term
        Terminal aminoacids that should be kept.
    """
    def __init__(
        self,
        cut: str,
        nocut: str | None = None,
        sense: Literal['N', 'C'] = 'C',
        keep_term: Literal['N', 'C', 'both', None] = None,
    ) -> None:
        # A lot of type-guarding...
        if not isinstance(cut, str):
            raise TypeError("Cut aminoacids must be string")
        if not cut:
            raise ValueError("Need string for cut aminoacids")
        for aa in cut:
            if aa not in AMINOACIDS:
                raise ValueError(
                    f"Not an standard aminoacid single-letter code: '{aa}'"
                )

        if not isinstance(nocut, str | None):
            raise TypeError("No-cut aminoacids must be string or None")
        if nocut is not None:
            if not nocut:
                raise ValueError("Need string no-cut aminoacids (or None)")
            for aa in nocut:
                if aa not in AMINOACIDS:
                    raise ValueError(
                        f"Not an standard aminoacid single-letter code: '{aa}'"
                    )

        if nocut is not None and (shared := set(cut) & set(nocut)):
            raise ValueError(f"Shared cut and nocut aminoacids: {"".join(shared)}")

        if not isinstance(sense, str) or not sense or sense not in 'NC':
            raise TypeError("Cleavage sense must be 'N' or 'C'")

        match keep_term:
            case 'N':
                pattern = rf"(^.|[{cut}])"
            case 'C':
                pattern = rf"([{cut}]|.$)"
            case 'both':
                pattern = rf"(^.|[{cut}]|.$)"
            case None:
                pattern = rf"([{cut}])"
            case _:
                raise TypeError(
                    f"Kept termini must be 'N', 'C', 'both' or None, not '{keep_term}'"
                )

        if nocut is not None:
            pattern += rf"(?!{nocut})"

        self.__cut = cut
        self.__nocut = nocut
        self.__sense: Literal['N', 'C'] = sense
        self.__keep_term: Literal['N', 'C', 'both', None] = keep_term

        self._pattern = re.compile(pattern)

    @overload
    def __call__(self, sequence: Seq) -> Seq: ...

    @overload
    def __call__(self, sequence: MutableSeq) -> MutableSeq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

    @abstractmethod
    def __call__(self, sequence: SeqLike) -> SeqLike:
        """Receive a sequence and return a decoy based on enzymatic peptides.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        A version of `sequence`, according to the enzyme specifications given
        at class instantiation.
        """
        pass

    @property
    def cut(self) -> str:
        """Cleavage sites as a string."""
        return self.__cut

    @property
    def sense(self) -> Literal['N', 'C']:
        """Sense of cleavage."""
        return self.__sense

    @property
    def nocut(self) -> str | None:
        """Aminoacids that stop cleavage as a string."""
        return self.__nocut

    @property
    def keep_term(self) -> Literal['N', 'C', 'both', None]:
        """Terminal aminoacids that should be kept."""
        return self.__keep_term


class ReversePep(EnzymeSpecificGenerator):
    """Appliy pseudo-reverse decoy generation with the specified enzyme
    properties.

    Pseudo-reverse (or reverse peptide) means that the enzymatic peptides will
    be reversed, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> TPKYSQRQHT`

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

    Parameters
    ----------
    cut
        Cleavage sites as a string.
    sense
        Whether the enzyme cleaves the 'C' or 'N' bond of the cleavage site.
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids followed by these.
    keep_term
        Terminal aminoacids that should be kept.

    Examples
    --------
    >>> rev = ReversePep("R", sense="N")
    >>> print(rev.cut, rev.nocut, rev.sense, rev.keep_term, sep=', ')
    R, None, N, None
    >>> rev = ReversePep("KR", nocut="P", keep_term="both")
    >>> print(rev.cut, rev.nocut, rev.sense, rev.keep_term, sep=', ')
    KR, P, C, both

    Cut argument cannot be an empty string:

    >>> rev = ReversePep("")
    Traceback (most recent call last):
        ...
    ValueError: Need string for cut aminoacids

    Aminoacids must be one of the 20 standard aminoacid single-letter codes:

    >>> rev = ReversePep("KR", nocut="B")
    Traceback (most recent call last):
        ...
    ValueError: Not an standard aminoacid single-letter code: 'B'
    """

    @overload
    def __call__(self, sequence: Seq) -> Seq: ...

    @overload
    def __call__(self, sequence: MutableSeq) -> MutableSeq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

    @override
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
        >>> rev = ReversePep("KR", nocut="P")
        >>> rev('QSYKPTRTHQ')
        'TPKYSQRQHT'
        >>> rev = ReversePep("K", sense="N", keep_term="N")
        >>> rev('QSYKPTRTHQ')
        'QYSKQHTRTP'
        """
        fragments = re.split(self._pattern, str(sequence))
        rev_frags = [frag[::-1] for frag in fragments]
        return seq_cast(sequence, "".join(rev_frags))

    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to ``ReversePep(sequence)`` where
        `sequence` is a `str`.
        """
        fragments = re.split(self._pattern, str(sequence))
        rev_frags = [frag[::-1] for frag in fragments]
        return "".join(rev_frags)

    @overload
    def decoy_from_Seq(self, sequence: Seq) -> Seq: ...

    @overload
    def decoy_from_Seq(self, sequence: MutableSeq) -> MutableSeq: ...

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq | MutableSeq:
        """Convenience funcion. Equivalent to ``ReversePep(sequence)`` where
        `sequence` is a `Seq` or `MutableSeq`.
        """
        from .. import _bio
        _bio.register()
        return self.decoy_from_Seq(sequence)


class ShufflePep(EnzymeSpecificGenerator):
    """Appliy pseudo-shuffle decoy generation with the specified enzyme
    properties.

    Pseudo-shuffle (or shuffle peptide) means that the enzymatic peptides will
    be shuffled, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> YTSKQPRQHT`

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

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
    >>> shuf = ShufflePep("R", sense="N")
    >>> print(shuf.cut, shuf.nocut, shuf.sense, shuf.keep_term, sep=', ')
    R, None, N, None
    >>> shuf = ShufflePep("KR", nocut="P", keep_term="both")
    >>> print(shuf.cut, shuf.nocut, shuf.sense, shuf.keep_term, sep=', ')
    KR, P, C, both

    Cut argument cannot be an empty string:

    >>> shuf = ShufflePep("")
    Traceback (most recent call last):
        ...
    ValueError: Need string for cut aminoacids

    Aminoacids must be one of the 20 standard aminoacid single-letter codes:

    >>> shuf = ShufflePep("KR", nocut="B")
    Traceback (most recent call last):
        ...
    ValueError: Not an standard aminoacid single-letter code: 'B'
    """

    @overload
    def __call__(self, sequence: Seq) -> Seq: ...

    @overload
    def __call__(self, sequence: MutableSeq) -> MutableSeq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

    @override
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
        >>> shuf = ShufflePep("KR", nocut="P")
        >>> shuf('QSYKPTRTHQ')
        'YTSKQPRQHT'
        >>> shuf = ShufflePep("K", sense="N", keep_term="N")
        >>> shuf('QSYKPTRTHQ')
        'QSYKTHQPTR'
        """
        fragments = re.split(self._pattern, str(sequence))

        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return seq_cast(sequence, "".join(shuf_frags))

    def decoy_from_str(self, sequence: str) -> str:
        """Convenience funcion. Equivalent to ``ShufflePep(sequence)`` where
        `sequence` is a `str`.
        """
        fragments = re.split(self._pattern, str(sequence))
        shuf_frags = [self._shuffle(frag) for frag in fragments]
        return "".join(shuf_frags)

    @overload
    def decoy_from_Seq(self, sequence: Seq) -> Seq: ...

    @overload
    def decoy_from_Seq(self, sequence: MutableSeq) -> MutableSeq: ...

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq | MutableSeq:
        """Convenience funcion. Equivalent to ``ShufflePep(sequence)`` where
        `sequence` is a `Seq` or `MutableSeq`.
        """
        from .. import _bio
        _bio.register()
        return self.decoy_from_Seq(sequence)

    def _shuffle(self, frag: str) -> str:
        new = list(frag)
        RAND.shuffle(new)
        return "".join(new)


# Hackish solution, but it allows the code to always return the correct type
# without importing Biopython or deferring to another module
def seq_cast(obj: SeqLike, sequence: str) -> SeqLike:
    """Convenience function. Transforms a `sequence` str into the correct
    :obj:`SeqLike` representation (through `obj`).

    This function doesn't need Biopython installed.

    Examples
    --------
    >>> from Bio.Seq import Seq, MutableSeq
    >>> str_seq = 'QSYKPTRTHQ'
    >>> bio_seq = Seq('YTSKQPRQHT')
    >>> seq_cast(bio_seq, str_seq)
    Seq('QSYKPTRTHQ')
    >>> bio_seq = MutableSeq('YTSKQPRQHT')
    >>> seq_cast(bio_seq, str_seq)
    MutableSeq('QSYKPTRTHQ')
    """
    cls = type(obj)
    return cls(sequence)

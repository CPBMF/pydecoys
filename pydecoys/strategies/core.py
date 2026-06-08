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

import random
import re
import typing as t
from abc import ABC, abstractmethod

if t.TYPE_CHECKING:
    from Bio.Seq import MutableSeq, Seq


type SeqLike = 'str | Seq | MutableSeq'
"""`SeqLike` objects can be indexed and spliced; `str` at runtime."""

type DecoyGenerator[T: SeqLike] = t.Callable[[T], T]
"""TypeAlias specifying the signature for decoy strategies.

A decoy strategy should be a ``Callable[[T], T]`` where ``T`` is a
:type:`SeqLike`.
"""

# So shuffled decoys are always reproducible
RAND: t.Final = random.Random(10)
"""Random number generator for stochastic decoy strategies."""

STD_AMINOACIDS: t.Final = 'QWERTYIPASDFGHKLCVNM'
"""Standard 20 aminoacids single-letter codes, majuscule."""

EXT_AMINOACIDS: t.Final = STD_AMINOACIDS + 'OU' + 'BJZX'
"""Extended aminoacids single-letter codes, majuscule.

Non-standard aminoacids
-----------------------

- O: Pyrrolysine
- U: Selenocysteine

Special codes
-------------

- B: Aspartic acid (D) or Asparagine (N)
- J: Leucine (L) or Isoleucine (I)
- Z: Glutamic acid (E) or Glutamine (Q)
- X: Any aminoacid

Notes
-----

Special letter codes are treated as literal characters when matching against
proteins. For example, 'B' won't match against either 'D' or 'N'. Those are
meant as ambiguous stand-ins in fasta protein sequences, not as regex
shortcuts.
"""


@t.runtime_checkable
class ContextfulGenerator(t.Protocol):
    """Protocol defining a decoy generator function that uses previously
    learned context.
    """

    def learn_context(self, sequences: t.Iterable[SeqLike]) -> None:
        """Receive the target proteins set to generate the necessary context.

        Parameters
        ----------
        sequences
            The target dataset.
        """
        ...

    def reset(self) -> None:
        """Reset the generator, so it holds no context and :attr:`is_set` i
        `False`.
        """
        ...

    @property
    def is_set(self) -> bool:
        """Whether the generator has context (`True`) or not (`False`)."""
        ...

    def __call__[T: SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return a decoy based on previous context.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A decoy version of `sequence`.
        """
        ...


class EnzymeSpecificGenerator(ABC):
    """Abstract base class for enzymatic aware decoy generation.

    Sequences can be split with :meth:`split_sequence`. The class is
    case-insensitive.

    The regex MUST match **only** the cleavage sites that shouldn't be
    altered. The cleavage sites MUST be captured by the regex pattern. Else,
    the resulting iterator from :meth:`split_sequence` won't yield all
    aminoacid residues.

    Parameters
    ----------
    pattern
        A regex pattern that must capture the desired cleavage sites. For
        example, for trypsin: ``r'([KR])(?!P)'``.
    sense
        Whether the enzyme cleaves the C-terminal, N-terminal or both termini
        of the cleavage site. This is unused by default, but can be useful for
        subclasses overriding the class. Case sensitive.

    Examples
    --------
    >>> class DummyEnzymeGenerator(EnzymeSpecificGenerator):
    ...     def __call__(sequence): raise NotImplementedError
    >>> dummy = DummyEnzymeGenerator(r"(R)", sense="N")
    >>> print(dummy.pattern)
    re.compile('(R)', re.IGNORECASE)
    >>> dummy = DummyEnzymeGenerator(r"([KR])(?!P)")
    >>> print(dummy.pattern)
    re.compile('([KR])(?!P)', re.IGNORECASE)

    The pattern argument cannot be an empty string:

    >>> dummy = DummyEnzymeGenerator("")
    Traceback (most recent call last):
        ...
    ValueError: Need string or re.Pattern for pattern
    """

    def __init__(
        self,
        pattern: str | re.Pattern[str],
        sense: t.Literal['N', 'C', 'both'] = 'C'
    ):
        if isinstance(pattern, str):
            if not pattern:
                raise ValueError('Need string or re.Pattern for pattern')
            self.__pattern = re.compile(pattern, re.IGNORECASE)
        elif isinstance(pattern, re.Pattern):
            self.__pattern = pattern
        else:
            raise TypeError('The pattern must be a string or a re.Pattern obj')

    def split_sequence(
        self,
        sequence: SeqLike
    ) -> t.Generator[tuple[str, bool], None, None]:
        """Split a given sequence into enzymatic fragments (minus the clevage
        site) and cleavage sites, in the order they appear.

        Parameters
        ----------
        sequence
            Aminoacid sequence to be split.

        Yields
        ------
        tuple[str, bool]
            A tuple containin an enzymatic fragment (minus the clevage site)
            and `False`, or a cleavage site and `True`.

        Examples
        --------
        >>> class DummySpecificGenerator(EnzymeSpecificGenerator):
        ...     def __call__(self, sequence): raise NotImplementedError
        >>> splitter = DummySpecificGenerator.from_enzyme('KR')
        >>> for val in splitter.split_sequence('QSYKPTRTHQ'):
        ...     print(val)
        ('QSY', False)
        ('K', True)
        ('PT', False)
        ('R', True)
        ('THQ', False)
        """

        for i, frag in enumerate(re.split(self.__pattern, str(sequence))):
            if frag:
                # Captured values (in this case, cleavage sites) are
                # guaranteed to be in odd indexes
                yield frag, i % 2 == 1

    @abstractmethod
    def __call__[T: SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return a decoy based on enzymatic peptides.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A decoy version of `sequence`, according to the enzyme
            specifications given at class instantiation.
        """
        ...

    @property
    def pattern(self) -> re.Pattern:
        """Regex pattern to capture cleavage sites."""
        return self.__pattern

    @classmethod
    def from_enzyme(
        cls,
        cut: str,
        nocut: str | None = None,
        nocut_n: str | None = None,
        sense: t.Literal['N', 'C', 'both'] = 'C',
    ) -> t.Self:
        """
        Create regex pattern from enzyme specifications and return a new
        instance with this regex pattern.

        Parameters
        ----------
        cut
            Cleavage sites as a string.
        nocut
            Aminoacids that stop cleavage when at C-terminal as a string, or
            `None`. If given, the enzyme will ignore `cut` aminoacids followed
            by these at the C-terminal.
        nocut_n
            Aminoacids that stop cleavage when at N-terminal as a string, or
            `None`. If given, the enzyme will ignore `cut` aminoacids
            preceeded by these at the N-terminal.
        sense
            Whether the enzyme cleaves the C-terminal, N-terminal or both
            termini of the cleavage site. This is unused by default, but can
            be useful for subclasses overriding the class. Case sensitive.

        Returns
        -------
        Self
            An instance of the class with a regex pattern constructed from the
            parameters as such: ``rf'(?<![{nocut_n}])([{cut}])(?![{nocut}])'``.

        Examples
        --------
        >>> class DummyEnzymeGenerator(EnzymeSpecificGenerator):
        ...     def __call__(sequence): raise NotImplementedError
        >>> dummy = DummyEnzymeGenerator.from_enzyme("R", sense="N")
        >>> print(dummy.pattern)
        re.compile('(R)', re.IGNORECASE)
        >>> dummy = DummyEnzymeGenerator.from_enzyme("KR", nocut="P")
        >>> print(dummy.pattern)
        re.compile('([KR])(?!P)', re.IGNORECASE)

        Cut argument cannot be an empty string:

        >>> dummy = DummyEnzymeGenerator.from_enzyme("")
        Traceback (most recent call last):
            ...
        ValueError: Need string for cut aminoacids

        Aminoacids must be one of the :data:`EXT_AMINOACIDS` single-letter
        codes:

        >>> dummy = DummyEnzymeGenerator.from_enzyme("KR", nocut="7")
        Traceback (most recent call last):
            ...
        ValueError: Not a valid aminoacid single-letter code: '7'
        """

        # A lot of type-guarding...
        if not isinstance(cut, str):
            raise TypeError("Cut aminoacids must be string")
        cut = cut.upper()
        if not cut:
            raise ValueError("Need string for cut aminoacids")
        cls._check_if_aa(cut)

        if nocut is None:
            pass
        elif isinstance(nocut, str):
            if not nocut:
                raise ValueError("Need string for nocut aminoacids (or None)")
            nocut = nocut.upper()
            cls._check_if_aa(nocut)
        else:
            raise TypeError("Nocut aminoacids must be string or None")

        if nocut_n is None:
            pass
        elif isinstance(nocut_n, str):
            if not nocut_n:
                raise ValueError("Need string for nocut_n aminoacids (or None)")
            nocut_n = nocut_n.upper()
            cls._check_if_aa(nocut_n)
        else:
            raise TypeError("Nocut_n aminoacids must be string or None")

        pattern = rf"([{cut}])" if len(cut) > 1 else rf"({cut})"

        if nocut is not None:
            pattern += rf"(?![{nocut}])" if len(nocut) > 1 else rf"(?!{nocut})"
        if nocut_n is not None:
            pattern = (
                (
                    rf"(?<![{nocut_n}])"
                    if len(nocut_n) > 1
                    else rf"(?<!{nocut_n})"
                )
                + pattern
            )

        return cls(pattern, sense)

    @staticmethod
    def _check_if_aa(sequence: str):
        """Raise a ValueError if a given sequence is not composed of
        :data:`EXT_AMINOACIDS` only.
        """

        for aa in sequence:
            if aa not in EXT_AMINOACIDS:
                raise ValueError(f"Not a valid aminoacid single-letter code: '{aa}'")


# Hackish solution, but it allows the code to always return the correct type
# without importing Biopython or deferring to another module
def seq_cast[T: SeqLike](obj: T, sequence: str) -> T:
    """Convenience function. Transforms a `sequence` str into the correct
    :type:`SeqLike` representation (through `obj`).

    This function doesn't need Biopython installed.

    Parameters
    ----------
    obj
        A ``SeqLike`` object specifying the return type.
    sequence
        A str sequence.

    Returns
    -------
    T
        The `sequence` as the type of `obj`.

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
    return cls(sequence)  # type: ignore

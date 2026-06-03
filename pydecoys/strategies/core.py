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
        A tuple containin an enzymatic fragment (minus the clevage site) and
        `False`, or a cleavage site and `True`.

        Examples
        --------
        >>> class DummySpecificGenerator(EnzymeSpecificGenerator):
        ...     def __call__(self, sequence): raise NotImplementedError
        >>> splitter = DummySpecificGenerator('KR')
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
        A decoy version of `sequence`, according to the enzyme specifications
        given at class instantiation.
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

        Examples
        --------
        >>> class DummyEnzymeGenerator(EnzymeSpecificGenerator):
        ...     def __call__(sequence): raise NotImplementedError
        >>> dummy = DummyEnzymeGenerator("R", sense="N")
        >>> print(dummy.pattern)
        re.compile('([R])', re.IGNORECASE)
        >>> dummy = DummyEnzymeGenerator("KR", nocut="P")
        >>> print(dummy.pattern)
        re.compile('([KR])(?![P])', re.IGNORECASE)

        Cut argument cannot be an empty string:

        >>> dummy = DummyEnzymeGenerator("")
        Traceback (most recent call last):
            ...
        ValueError: Need string for cut aminoacids

        Aminoacids must be one of the :data:`EXT_AMINOACIDS` single-letter
        codes:

        >>> dummy = DummyEnzymeGenerator("KR", nocut="7")
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

        pattern = rf"([{cut}])"

        if nocut is not None:
            pattern += rf"(?![{nocut}])"
        if nocut_n is not None:
            pattern = rf"(?<![{nocut_n}])" + pattern

        return cls(pattern, sense)

    @staticmethod
    def _check_if_aa(sequence: str):
        """Raise a ValueError if a given sequence is not composed of
        :data:`EXT_AMINOACIDS` only.
        """

        for aa in sequence:
            if aa not in EXT_AMINOACIDS:
                raise ValueError(f"Not a valid aminoacid single-letter code: '{aa}'")


class ReversePep(EnzymeSpecificGenerator):
    """Apply pseudo-reverse decoy generation with the specified enzyme
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
    nocut
        Aminoacids that stop cleavage at C-bond as a string, or `None`. If
        given, the enzyme will ignore `cut` aminoacids followed by these at
        the C position.
    nocut_n
        Aminoacids that stop cleavage at N-bond as a string, or `None`. If
        given, the enzyme will ignore `cut` aminoacids preceeded by these at
        the N position.
    sense
        Whether the enzyme cleaves the C or N bond of the cleavage site.
    """

    @t.override
    def __call__[T: SeqLike](self, sequence: T) -> T:
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
        >>> rev = ReversePep("K", sense="N")
        >>> rev('QSYKPTRTHQ')
        'YSQKQHTRTP'
        """

        fragments = self.split_sequence(sequence)
        rev_frags = []

        for frag, cleavage in fragments:
            if not cleavage:
                frag = frag[::-1]
            rev_frags.append(frag)

        return seq_cast(sequence, "".join(rev_frags))


class ShufflePep(EnzymeSpecificGenerator):
    """Apply pseudo-shuffle decoy generation with the specified enzyme
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
    nocut
        Aminoacids that stop cleavage at C-bond as a string, or `None`. If
        given, the enzyme will ignore `cut` aminoacids followed by these at
        the C position.
    nocut_n
        Aminoacids that stop cleavage at N-bond as a string, or `None`. If
        given, the enzyme will ignore `cut` aminoacids preceeded by these at
        the N position.
    sense
        Whether the enzyme cleaves the C or N bond of the cleavage site.
    """

    @t.override
    def __call__[T: SeqLike](self, sequence: T) -> T:
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
        >>> shuf = ShufflePep("K", sense="N")
        >>> shuf('QSYKPTRTHQ')
        'QYSKTHRPTQ'
        """

        fragments = self.split_sequence(sequence)
        shuf_frags = []

        for frag, cleavage in fragments:
            if not cleavage:
                frag = self._shuffle(frag)
            shuf_frags.append(frag)

        return seq_cast(sequence, "".join(shuf_frags))

    @staticmethod
    def _shuffle(frag: str) -> str:
        new = list(frag)
        RAND.shuffle(new)
        return "".join(new)


class RandomizePep(EnzymeSpecificGenerator):
    """Apply pseudo-randomize decoy generation with the specified enzyme
    properties.

    Pseudo-randomize (or randomize peptide) means that the enzymatic peptides
    will be randomized, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> DSDPCCRGIS`

    The randomization happens based on aminoacid proportions from the target
    database. Cleavage sites aren't counted.

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

    Parameters
    ----------
    cut
        Cleavage sites as a string.
    nocut
        Aminoacids that stop cleavage at C-bond as a string, or `None`. If
        given, the enzyme will ignore `cut` aminoacids followed by these at
        the C position.
    nocut_n
        Aminoacids that stop cleavage at N-bond as a string, or `None`. If
        given, the enzyme will ignore `cut` aminoacids preceeded by these at
        the N position.
    sense
        Whether the enzyme cleaves the C or N bond of the cleavage site.
    """

    _AA_TO_INDEX = {aa: i for i, aa in enumerate(EXT_AMINOACIDS)}

    @t.override
    def __init__(
        self,
        pattern: str | re.Pattern[str],
        sense: t.Literal['N', 'C', 'both'] = 'C',
    ) -> None:
        super().__init__(pattern, sense)
        self._weights: list[int] | None = None

    @t.override
    def __call__[T: SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return a pseudo-randomized decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        A pseudo-randomized version of `sequence`, according to the enzyme
        specifications given at class instantiation.

        Examples
        --------
        >>> rand = RandomizePep("KR", nocut="P")
        >>> rand('QSYKPTRTHQ')  # doctest: +SKIP
        'DSDPCCRGIS'
        >>> rand = ShufflePep("K", sense="N")
        >>> rand('QSYKPTRTHQ')  # doctest: +SKIP
        'PINKMEVDAP'
        """

        rand_frags = []
        fragments = self.split_sequence(sequence)

        for frag, cleavage in fragments:
            if not cleavage:
                frag = self._get_rand(frag)
            rand_frags.append(frag)

        decoy = "".join(rand_frags)
        return seq_cast(sequence, decoy)

    def learn_context(self, sequences: t.Iterable[SeqLike]) -> None:
        """Receive the target proteins set to learn aminoacid proportions and
        use them as weights during randomization.

        Since cleavage sites are unaltered during randomization, they are
        ignored here, so proportions are kept equal.

        Any character not in :data:`EXT_AMINOACIDS` is ignored.

        Parameters
        ----------
        sequences
            The target dataset.
        """

        self._weights = [0] * len(EXT_AMINOACIDS)

        for seq in sequences:
            for frag, cleavage in self.split_sequence(seq):
                # We don't count cleavage sites in the weights since they'll
                # be directly preserved
                if cleavage:
                    continue
                for aa in frag:
                    idx = self._AA_TO_INDEX.get(aa.upper())
                    if idx is not None:
                        self._weights[idx] += 1

    def reset(self) -> None:
        """Reset the generator, erasing its previous context."""
        self._weights = None

    @property
    def is_set(self) -> bool:
        """Whether the generator has context (`True`) or not (`False`)."""
        return self._weights is not None

    def _get_rand(self, frag: str) -> str:
        length = len(frag)
        new = RAND.choices(EXT_AMINOACIDS, weights=self._weights, k=length)
        return "".join(new)


# Hackish solution, but it allows the code to always return the correct type
# without importing Biopython or deferring to another module
def seq_cast[T: SeqLike](obj: T, sequence: str) -> T:
    """Convenience function. Transforms a `sequence` str into the correct
    :type:`SeqLike` representation (through `obj`).

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
    return cls(sequence)  # type: ignore

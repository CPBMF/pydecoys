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

    This class creates a compiled regex pattern at instantiation that captures
    cleavage sites. A sequence can be split with :meth:`split_sequence`.

    This class also save the enzymatic specifications as get-only attributes.

    The nocut value is always considered at the C-terminal, even if the sense
    is N.

    Parameters
    ----------
    cut
        Cleavage sites as a string.
    sense
        Whether the enzyme cleaves the C or N bond of the cleavage site.
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids with a C-terminal followed by these.
        The nocut value is always at the C-terminal, even if the sense is N.
    """
    def __init__(
        self,
        cut: str,
        nocut: str | None = None,
        sense: t.Literal['N', 'C'] = 'C',
    ) -> None:
        # A lot of type-guarding...
        if not isinstance(cut, str):
            raise TypeError("Cut aminoacids must be string")
        if not cut:
            raise ValueError("Need string for cut aminoacids")
        for aa in cut:
            if aa not in EXT_AMINOACIDS:
                raise ValueError(f"Not a valid aminoacid single-letter code: '{aa}'")

        if not isinstance(nocut, str | None):
            raise TypeError("No-cut aminoacids must be string or None")
        if nocut is not None:
            if not nocut:
                raise ValueError("Need string no-cut aminoacids (or None)")
            for aa in nocut:
                if aa not in EXT_AMINOACIDS:
                    raise ValueError(
                        f"Not a valid aminoacid single-letter code: '{aa}'"
                    )

        if nocut is not None and (shared := set(cut) & set(nocut)):
            raise ValueError(f"Shared cut and nocut aminoacids: {"".join(shared)}")

        if not isinstance(sense, str) or not sense or sense not in {'N', 'C'}:
            raise TypeError("Cleavage sense must be 'N' or 'C'")

        pattern = rf"([{cut}])"

        if nocut is not None:
            pattern += rf"(?!{nocut})"

        self.__cut = cut
        self.__nocut = nocut
        self.__sense: t.Literal['N', 'C'] = sense
        self.__pattern = re.compile(pattern)

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
        A tuple containin an enzymatic fragments (minus the clevage site) and
        `False`, or a cleavage site and `True`. Cleavage sites are guaranteed
        to be one character only.

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
    def cut(self) -> str:
        """Cleavage sites as a string."""
        return self.__cut

    @property
    def sense(self) -> t.Literal['N', 'C']:
        """Sense of cleavage."""
        return self.__sense

    @property
    def nocut(self) -> str | None:
        """Aminoacids that stop cleavage as a string."""
        return self.__nocut


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
    sense
        Whether the enzyme cleaves the 'C' or 'N' bond of the cleavage site.
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids followed by these.

    Examples
    --------
    >>> rev = ReversePep("R", sense="N")
    >>> print(rev.cut, rev.nocut, rev.sense, sep=', ')
    R, None, N
    >>> rev = ReversePep("KR", nocut="P")
    >>> print(rev.cut, rev.nocut, rev.sense, sep=', ')
    KR, P, C

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

        # Cleavage sites are guaranteed to always be one letter only,
        # reverting them is no-op
        fragments = self.split_sequence(sequence)
        rev_frags = [frag[0][::-1] for frag in fragments]
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
    sense
        Sense cleavage (whether the enzyme cleaves the 'C' or 'N' bond of the
        cleavage site).
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids followed by these.

    Examples
    --------
    >>> shuf = ShufflePep("R", sense="N")
    >>> print(shuf.cut, shuf.nocut, shuf.sense, sep=', ')
    R, None, N
    >>> shuf = ShufflePep("KR", nocut="P")
    >>> print(shuf.cut, shuf.nocut, shuf.sense, sep=', ')
    KR, P, C

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

        # Cleavage sites are guaranteed to always be one letter only,
        # shuffling them is no-op
        fragments = self.split_sequence(sequence)
        shuf_frags = [self._shuffle(frag[0]) for frag in fragments]
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
    sense
        Whether the enzyme cleaves the 'C' or 'N' bond of the cleavage site.
    nocut
        Aminoacids that stop cleavage as a string, or `None`. If given, the
        enzyme won't cut aminoacids followed by these.

    Examples
    --------
    >>> rand = RandomizePep("R", sense="N")
    >>> print(rand.cut, rand.nocut, rand.sense, sep=', ')
    R, None, N
    >>> rand = RandomizePep("KR", nocut="P")
    >>> print(rand.cut, rand.nocut, rand.sense, sep=', ')
    KR, P, C

    Cut argument cannot be an empty string:

    >>> rand = RandomizePep("")
    Traceback (most recent call last):
        ...
    ValueError: Need string for cut aminoacids

    Aminoacids must be one of the 20 standard aminoacid single-letter codes:

    >>> rand = RandomizePep("KR", nocut="B")
    Traceback (most recent call last):
        ...
    ValueError: Not an standard aminoacid single-letter code: 'B'
    """

    _AA_TO_INDEX = {aa: i for i, aa in enumerate(EXT_AMINOACIDS)}

    @t.override
    def __init__(
        self,
        cut: str,
        nocut: str | None = None,
        sense: t.Literal['N', 'C'] = 'C',
    ) -> None:
        super().__init__(cut, nocut, sense)
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
                    idx = self._AA_TO_INDEX.get(aa)
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

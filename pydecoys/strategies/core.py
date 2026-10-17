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

# TODO: Fix docstrings

"""Core API (protocols, ABCs, types and utility) for :mod:`strategies`."""

from __future__ import annotations

import random
import re
import typing as t
from abc import ABC, abstractmethod
from collections import Counter, defaultdict

if t.TYPE_CHECKING:
    from Bio.Seq import MutableSeq, Seq


type SeqLike = 'str | Seq | MutableSeq'
"""`SeqLike` objects can be indexed and spliced; `str` at runtime."""

type SimpleGenerator[T: SeqLike] = t.Callable[[T], T]
"""TypeAlias specifying the signature for decoy strategies that don't take
cleavage agents into account.
"""

type CleavageGenerator[T: SeqLike] = t.Callable[
    [T, re.Pattern, t.Literal['N', 'C', 'both']],
    T
]
"""TypeAlias specifying the signature for decoy strategies that take cleavage
agents into account.
"""

type DecoyGenerator = SimpleGenerator | CleavageGenerator
"""TypeAlias specifying the signature for all decoy strategies.
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

    def __call__[T: SeqLike](self, sequence: T, *args) -> T:
        """Receive a sequence and return a decoy based on previous context.

        Parameters
        ----------
        sequence
            A single sequence.
        args
            Arguments for cleavage-aware generators.

        Returns
        -------
        T
            A decoy version of `sequence`.
        """
        ...


class EnzymeSpecificGenerator(ABC):
    """Abstract base class for enzymatic aware decoy generation.

    Sequences can be split with :meth:`split_sequence`.
    """
    def split_sequence(
        self,
        sequence: SeqLike,
        pattern: re.Pattern
    ) -> t.Generator[tuple[str, bool], None, None]:
        """Split a given sequence into enzymatic fragments (minus the clevage
        site) and cleavage sites, in the order they appear.

        The regex MUST match **only** the cleavage sites that shouldn't be
        altered. The cleavage sites MUST be captured by the regex pattern. Else,
        the resulting iterator from :meth:`split_sequence` won't yield all
        aminoacid residues.

        The function is case-insensitive.

        Parameters
        ----------
        sequence
            Aminoacid sequence to be split.
        pattern
            A regex pattern that matches cleavage sites.

        Yields
        ------
        tuple[str, bool]
            A tuple containin an enzymatic fragment (minus the clevage site)
            and `False`, or a cleavage site and `True`.

        Examples
        --------
        >>> class DummySpecificGenerator(EnzymeSpecificGenerator):
        ...     def __call__(self, sequence): raise NotImplementedError
        >>> splitter = DummySpecificGenerator()
        >>> for val in splitter.split_sequence('QSYKPTRTHQ', r'([KR])'):
        ...     print(val)
        ('QSY', False)
        ('K', True)
        ('PT', False)
        ('R', True)
        ('THQ', False)
        """

        for i, frag in enumerate(re.split(pattern, str(sequence), flags=re.IGNORECASE)):
            if frag:
                # Captured values (in this case, cleavage sites) are
                # guaranteed to be in odd indexes
                yield frag, i % 2 == 1

    @abstractmethod
    def __call__[T: SeqLike](
        self,
        sequence: T,
        pattern: re.Pattern,
        sense = t.Literal['N', 'C', 'both']
    ) -> T:
        """Receive a sequence and return a decoy based on enzymatic peptides.
        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A decoy version of `sequence`, according to the cleavage agent
            specification.
        """
        ...


class ReversePep(EnzymeSpecificGenerator):
    """Apply pseudo-reverse decoy generation with the specified enzyme
    properties.

    Pseudo-reverse (or reverse peptide) means that the enzymatic peptides will
    be reversed, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> TPKYSQRQHT`

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.
    """

    @t.override
    def __call__[T: SeqLike](
        self,
        sequence: T,
        pattern: re.Pattern,
        sense = t.Literal['N', 'C', 'both']
    ) -> T:
        fragments = self.split_sequence(sequence, pattern)
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
    """

    @t.override
    def __call__[T: SeqLike](
        self,
        sequence: T,
        pattern: re.Pattern,
        sense = t.Literal['N', 'C', 'both']
    ) -> T:
        fragments = self.split_sequence(sequence, pattern)
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
    """

    _AA_TO_INDEX = {aa: i for i, aa in enumerate(EXT_AMINOACIDS)}

    @t.override
    def __init__(self) -> None:
        self._weights: list[int] | None = None

    @t.override
    def __call__[T: SeqLike](
        self,
        sequence: T,
        pattern: re.Pattern,
        sense = t.Literal['N', 'C', 'both']
    ) -> T:
        rand_frags = []
        fragments = self.split_sequence(sequence, pattern)

        for frag, cleavage in fragments:
            if not cleavage:
                frag = self._get_rand(frag)
            rand_frags.append(frag)

        decoy = "".join(rand_frags)
        return seq_cast(sequence, decoy)

    def learn_context(
            self,
            sequences: t.Iterable[SeqLike],
            pattern: re.Pattern,
            sense: t.Literal['N', 'C', 'both']
        ) -> None:
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
            for frag, cleavage in self.split_sequence(seq, pattern):
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


class MarkovPep(EnzymeSpecificGenerator):
    """Apply pseudo-markov-chain decoy generation with the specified enzyme
    properties.

    Pseudo-markov-chain means that the enzymatic peptides will be randomized
    using a markov-chain model, except for the cleavage site.

    Cleavage sites aren't counted when computing the probability of going from
    a previous state to another, but are counted as previous state.

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

    Parameters
    ----------
    k
        The state-space for the Markov-chain model. The next aminoacid of the
        sequence will have its probability determined by the `k` aminoacids
        before.
    """

    @t.override
    def __init__(self, k: int) -> None:
        self._weights: dict[
            tuple[str | None, ...],
            tuple[list[str], list[int]]
        ] | None = None
        self._global_weights: tuple[list[str], list[int]] | None = None
        self._initial_state = (None, ) * k

    @t.override
    def __call__[T: SeqLike](
        self,
        sequence: T,
        pattern: re.Pattern,
        sense = t.Literal['N', 'C', 'both']
    ) -> T:
        """Receive a sequence and return a pseudo-markov-chain generated decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A pseudo-markov-chain generated version of `sequence`, according
            to the enzyme specifications given at class instantiation.

        Examples
        --------
        >>> rand = RandomizePep.from_enzyme("KR", nocut="P")
        >>> rand('QSYKPTRTHQ')  # doctest: +SKIP
        'DSDPCCRGIS'
        >>> rand = RandomizePep.from_enzyme("K", sense="N")
        >>> rand('QSYKPTRTHQ')  # doctest: +SKIP
        'PINKMEVDAP'
        """

        if self._weights is None:
            raise RuntimeError("The generator has no context.")

        markov_frags = []
        fragments = self.split_sequence(sequence, pattern)

        state: tuple[str | None, ...] = self._initial_state
        for frag, cleavage in fragments:
            if cleavage:
                markov_frags.append(frag)
                state = (*state[len(frag):], *frag)
                continue
            for aa in frag:
                try:
                    new = RAND.choices(
                        self._weights[state][0],
                        self._weights[state][1],
                        k=1
                    )[0]
                except IndexError:
                    new = RAND.choices(
                        self._global_weights[0],
                        self._global_weights[1],
                        k=1
                    )[0]
                markov_frags.append(new)
                state = (*state[1:], aa)

        decoy = "".join(markov_frags)
        return seq_cast(sequence, decoy)

    def learn_context(
        self,
        sequences: t.Iterable[SeqLike],
        pattern: re.Pattern,
        sense = t.Literal['N', 'C', 'both']
    ):
        weights: dict[tuple[str | None, ...], Counter[str]] = defaultdict(Counter)
        global_weights: Counter[str] = Counter()

        for sequence in sequences:
            prev: tuple[str | None, ...] = self._initial_state

            fragments = self.split_sequence(sequence, pattern)
            for frag, cleavage in fragments:
                if cleavage:
                    prev = (*prev[len(frag):], *frag)
                    continue
                for aa in frag:
                    weights[prev][aa] += 1
                    global_weights[aa] += 1
                    prev = (*prev[1:], aa)

        self._weights = {
            state: (list(counter.keys()), list(counter.values()))
            for state, counter
            in weights.items()
        }
        self._global_weights = (
            list(global_weights.keys()),
            list(global_weights.values())
        )

    def reset(self) -> None:
        """Reset the generator, erasing its previous context."""
        self._weights = None
        self._global_weights = None

    @property
    def is_set(self) -> bool:
        """Whether the generator has context (`True`) or not (`False`)."""
        return self._weights is not None


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

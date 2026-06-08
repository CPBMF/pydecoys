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

"""Internal logic for implemented strategies."""

import re
from collections import Counter, defaultdict
from typing import Iterable, Literal, cast, override

from pydecoys import _registry as r
from pydecoys import strategies as s
from pydecoys._registry import DECOY_STRATEGIES  # noqa: F401


@r.register_fn('reverse')
def reverse[T: s.SeqLike](sequence: T) -> T:
    """Return the reversed `sequence`, inverting C- and N-terminal positions.
    """
    return cast(T, sequence[::-1])


@r.register_fn('shuffle')
def shuffle[T: s.SeqLike](sequence: T) -> T:
    """Return the shuffled `sequence`, shuffling the aminoacids in place."""
    new = list(sequence)
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, "".join(new))


@r.register_cls('randomize')
class Randomize:
    _AA_TO_INDEX = {aa: i for i, aa in enumerate(s.EXT_AMINOACIDS)}

    def __init__(self):
        self._weights = None

    def learn_context(self, sequences: Iterable[s.SeqLike]):
        self._weights = [0] * len(s.EXT_AMINOACIDS)

        for seq in sequences:
            for aa in seq:
                idx = self._AA_TO_INDEX.get(aa.upper())
                if idx is not None:
                    self._weights[idx] += 1

    def reset(self) -> None:
        self._weights = None

    @property
    def is_set(self) -> bool:
        return self._weights is not None

    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        length = len(sequence)
        new = s.RAND.choices(s.EXT_AMINOACIDS, weights=self._weights, k=length)
        return s.seq_cast(sequence, "".join(new))


class Markov:
    def __init__(self, k: int = 1):
        self._weights: dict[
            tuple[str | None, ...],
            tuple[list[str], list[int]]
        ] | None = None
        self._global_weights: tuple[list[str], list[int]] = ([], [])
        self._initial_state = (None, ) * k

    def learn_context(self, sequences: Iterable[s.SeqLike]):
        weights: dict[tuple[str | None, ...], Counter[str]] = defaultdict(Counter)
        global_weights: Counter[str] = Counter()

        for sequence in sequences:
            prev: tuple[str | None, ...] = self._initial_state
            for aa in sequence:
                aa = aa.upper()
                if aa not in s.EXT_AMINOACIDS:
                    continue
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

    def reset(self):
        self._weights = None

    @property
    def is_set(self) -> bool:
        return self._weights is not None

    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        if self._weights is None:
            raise RuntimeError("The generator has no context.")

        prev: tuple[str | None, ...] = self._initial_state
        decoy = []
        for _ in sequence:
            try:
                aa = s.RAND.choices(
                    self._weights[prev][0],
                    self._weights[prev][1],
                    k=1
                )[0]
            except IndexError:
                # If we reach a state that has no following state, we'll
                # fallback to global distribution
                aa = s.RAND.choices(
                    self._global_weights[0],
                    self._global_weights[1],
                    k=1
                )[0]

            decoy.append(aa)
            prev = (*prev[1:], aa)

        return s.seq_cast(sequence, "".join(decoy))


# for i in range(1, 7):
#     r.add_callable(f'markov{i}', Markov(i))


@r.register_cleavage_aware('reversepep')
class ReversePep(s.EnzymeSpecificGenerator):
    """Apply pseudo-reverse decoy generation with the specified enzyme
    properties.

    Pseudo-reverse (or reverse peptide) means that the enzymatic peptides will
    be reversed, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> TPKYSQRQHT`

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

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
    """

    @override
    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return a pseudo-reversed decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A pseudo-reversed version of `sequence`, according to the enzyme
            specifications given at class instantiation.
        """

        fragments = self.split_sequence(sequence)
        rev_frags = []

        for frag, cleavage in fragments:
            if not cleavage:
                frag = frag[::-1]
            rev_frags.append(frag)

        return s.seq_cast(sequence, "".join(rev_frags))


@r.register_cleavage_aware('shufflepep')
class ShufflePep(s.EnzymeSpecificGenerator):
    """Apply pseudo-shuffle decoy generation with the specified enzyme
    properties.

    Pseudo-shuffle (or shuffle peptide) means that the enzymatic peptides will
    be shuffled, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> YTSKQPRQHT`

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

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
    """

    @override
    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return a pseudo-shuffled decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A pseudo-shuffled version of `sequence`, according to the enzyme
            specifications given at class instantiation.
        """

        fragments = self.split_sequence(sequence)
        shuf_frags = []

        for frag, cleavage in fragments:
            if not cleavage:
                frag = self._shuffle(frag)
            shuf_frags.append(frag)

        return s.seq_cast(sequence, "".join(shuf_frags))

    @staticmethod
    def _shuffle(frag: str) -> str:
        new = list(frag)
        s.RAND.shuffle(new)
        return "".join(new)


@r.register_cleavage_aware('randomizepep')
class RandomizePep(s.EnzymeSpecificGenerator):
    """Apply pseudo-randomize decoy generation with the specified enzyme
    properties.

    Pseudo-randomize (or randomize peptide) means that the enzymatic peptides
    will be randomized, except for the cleavage site. For trypsin:

    - `QSYKPTRTHQ -> QSYKPTR.THQ -> DSDPCCRGIS`

    Each aminoacid is given a weighted probability of being chosen based on
    its proportions from the target database. Cleavage sites aren't counted.

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

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
    """

    _AA_TO_INDEX = {aa: i for i, aa in enumerate(s.EXT_AMINOACIDS)}

    @override
    def __init__(
        self,
        pattern: str | re.Pattern[str],
        sense: Literal['N', 'C', 'both'] = 'C',
    ) -> None:
        super().__init__(pattern, sense)
        self._weights: list[int] | None = None

    @override
    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return a pseudo-randomized decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A pseudo-randomized version of `sequence`, according to the enzyme
            specifications given at class instantiation.
        """

        rand_frags = []
        fragments = self.split_sequence(sequence)

        for frag, cleavage in fragments:
            if not cleavage:
                frag = self._get_rand(frag)
            rand_frags.append(frag)

        decoy = "".join(rand_frags)
        return s.seq_cast(sequence, decoy)

    def learn_context(self, sequences: Iterable[s.SeqLike]) -> None:
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

        self._weights = [0] * len(s.EXT_AMINOACIDS)

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
        new = s.RAND.choices(s.EXT_AMINOACIDS, weights=self._weights, k=length)
        return "".join(new)


class MarkovPep(s.EnzymeSpecificGenerator):
    """Apply markov-chain decoy generation to enzymatic fragments, except
    cleavage sites.

    Cleavage sites aren't counted when computing the probability of a next
    state, but are counted as previous state.

    This better preserves actual peptide amount and sizes from the targets to
    the decoys.

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
    k
        The state-space for the Markov-chain model. The next aminoacid of the
        sequence will have its probability determined by the `k` aminoacids
        before it.
    """

    @override
    def __init__(
        self,
        pattern: str | re.Pattern[str],
        sense: Literal['N', 'C', 'both'] = 'C',
        k: int = 1,
    ) -> None:
        super().__init__(pattern, sense)
        self._weights: dict[
            tuple[str | None, ...],
            tuple[list[str], list[int]]
        ] | None = None
        self._initial_state = (None, ) * k

    @override
    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        """Receive a sequence and return an enzymatic-aware markov-chain
        generated decoy.

        Parameters
        ----------
        sequence
            A single sequence.

        Returns
        -------
        T
            A markov-chain generated version of `sequence`, according to the
            enzyme specifications given at class instantiation.
        """

        if self._weights is None:
            raise RuntimeError("The generator has no context.")

        markov_frags = []
        fragments = self.split_sequence(sequence)

        state: tuple[str | None, ...] = self._initial_state
        for frag, cleavage in fragments:
            if cleavage:
                markov_frags.append(frag)
                state = (*state[len(frag):], *frag)
                continue
            for aa in frag:
                try:
                    new = s.RAND.choices(
                        self._weights[state][0],
                        self._weights[state][1],
                        k=1
                    )[0]
                except IndexError:
                    new = s.RAND.choices(
                        self._global_weights[0],
                        self._global_weights[1],
                        k=1
                    )[0]
                markov_frags.append(new)
                state = (*state[1:], aa)

        decoy = "".join(markov_frags)
        return s.seq_cast(sequence, decoy)

    def learn_context(self, sequences: Iterable[s.SeqLike]):
        """Receive the target proteins set to learn aminoacid proportions and
        use them as weights during randomization.

        Since cleavage sites are unaltered during markov-chain generation,
        they are ignored as next state, but are still counted as previous
        state.

        Any character not in :data:`EXT_AMINOACIDS` is ignored.

        Parameters
        ----------
        sequences
            The target dataset.
        """

        weights: dict[tuple[str | None, ...], Counter[str]] = defaultdict(Counter)
        global_weights: Counter[str] = Counter()

        for sequence in sequences:
            prev: tuple[str | None, ...] = self._initial_state

            fragments = self.split_sequence(sequence)
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

    @property
    def is_set(self) -> bool:
        """Whether the generator has context (`True`) or not (`False`)."""
        return self._weights is not None


# Pre-defined enzymes
# Those specifications were taken from <https://github.com/HUPO-PSI/psi-ms-CV>
#
# (regex, cleavage sense, str key)
# The regex patterns MUST ONLY match the cleavage site and MUST capture it.
# Because of this, we cannot use the regex patterns from PSI MS Ontology.

type _Enzyme = tuple[str, Literal['N', 'C', 'both'], str]

_ENZYMES: list[_Enzyme] = [
    (r'([TASV])',        'C',    'alphalp'),        # AlphaLP
    (r'(R)(?!P)',        'C',    'argc'),           # Arg-C
    (r'([BD])',          'N',    'aspn'),           # Asp-N
    (r'([DE])',          'N',    'aspnambic'),      # Asp-N_ambic
    (r'([FYWL])(?!P)',   'C',    'chymo'),          # Chymotrypsin
    (r'(M)',             'C',    'cnbr'),           # CNBr
    (r'(D)',             'both', 'formicacid'),     # Formic_acid
    (r'(?<!E)(E)',       'C',    'gluc'),           # glutamyl endopeptidase
    (r'([ALIV])(?!P)',   'C',    'elastase'),       # leukocyte elastase
    (r'(K)(?!P)',        'C',    'lysc'),           # Lys-C
    (r'(K)',             'C',    'lyscp'),          # Lys-C/P
    (r'(K)',             'N',    'lysn'),           # Lys-N
    (r'([FL])',          'C',    'pepsina'),        # PepsinA
    (r'([HKR]P)(?!P)',   'C',    'proc'),           # proline endopeptidase
    (r'([KR])(?!P)',     'C',    'trypsin'),        # Trypsin
    (r'([KR])',          'C',    'trypsinp'),       # Trypsin/P
    (r'([FYWLKR])(?!P)', 'C',    'trypchymo'),      # TrypChymo
    (r'([KR])',          'N',    'trypn'),          # Tryp-N
    (r'(W)',             'C',    '2iodobenzoate'),  # 2-iodobenzoate
    (r'([BDEZ])(?!P)',   'C',    'v8de'),           # V8-DE
    (r'([EZ])(?!P)',     'C',    'v8e'),            # V8-E
]

for enzyme in _ENZYMES:
    r.add_cleavage_agent(enzyme[2], enzyme[0], enzyme[1])

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
from typing import Callable, Iterable, Literal, cast

from pydecoys import strategies as s


type _Enzyme = tuple[re.Pattern, Literal['N', 'C', 'both']]


decoy_strategies: dict[tuple[str, bool | str, bool], s.DecoyGenerator] = {}
"""The key is the strategy name, `True`, `False` or `'M'` on whether to keep
the N-terminal, and `True` or `False` on whether to keep the C-terminal.
"""


# Pre-defined enzymes
# Those specifications were taken from <https://github.com/HUPO-PSI/psi-ms-CV>
#
# (regex, cleavage sense)
# The regex patterns MUST ONLY match the cleavage site and MUST capture it.
# Because of this, we cannot use the regex patterns from PSI MS Ontology.

def _compile(pattern: str) -> re.Pattern:
    return re.compile(pattern)

enzymes: dict[str, _Enzyme] = {
    'AlphaLP': (_compile(r'([TASV])'), 'C'),
    'Arg-C': (_compile(r'(R)(?!P)'), 'C'),
    'Asp-N': (_compile(r'([BD])'), 'N'),
    'Asp-N_ambic': (_compile(r'([DE])'), 'N'),
    'Chymotrypsin': (_compile(r'([FYWL])(?!P)'), 'C'),
    'CNBr': (_compile(r'(M)'), 'C'),
    'Formic_acid': (_compile(r'(D)'), 'both'),
    'Glu-C': (_compile(r'(?<!E)(E)'), 'C'),
    'elastase': (_compile(r'([ALIV])(?!P)'), 'C'),
    'Lys-C': (_compile(r'(K)(?!P)'), 'C'),
    'Lys-C/P': (_compile(r'(K)'), 'C'),
    'Lys-N': (_compile(r'(K)'), 'N'),
    'PepsinA': (_compile(r'([FL])'), 'C'),
    'Pro-C': (_compile(r'([HKR]P)(?!P)'), 'C'),
    'Trypsin': (_compile(r'([KR])(?!P)'), 'C'),
    'Trypsin/P': (_compile(r'([KR])'), 'C'),
    'TrypChymo': (_compile(r'([FYWLKR])(?!P)'), 'C'),
    'Tryp-N': (_compile(r'([KR])'), 'N'),
    '2-iodobenzoate': (_compile(r'(W)'), 'C'),
    'V8-DE': (_compile(r'([BDEZ])(?!P)'), 'C'),
    'V8-E': (_compile(r'([EZ])(?!P)'), 'C'),
}


def reverse[T: s.SeqLike](sequence: T) -> T:
    """Return the reversed `sequence`."""
    return cast(T, sequence[::-1])


def shuffle[T: s.SeqLike](sequence: T) -> T:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, "".join(new))


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


_PLAIN_STRATEGIES: list[tuple[str, Callable[[], s.SimpleGenerator]]] = [
    ('reverse',   lambda: reverse),  # noqa: E272
    ('shuffle',   lambda: shuffle),  # noqa: E272
    ('randomize', Randomize)
]


_ENZYME_STRATEGIES: list[tuple[str, type[s.EnzymeSpecificGenerator]]] = [
    ('reversepep',   s.ReversePep),
    ('shufflepep',   s.ShufflePep),
    ('randomizepep', s.RandomizePep)
]


# Functional programming hell
def _register_all() -> None:
    for key, fn in _PLAIN_STRATEGIES:
        decoy_strategies[(key, False, False)] = fn()
        decoy_strategies[(key, True, False)] = s.keepsn(fn())
        decoy_strategies[(key, 'M', False)] = s.keepsm(fn())

        decoy_strategies[(key, False, True)] = s.keepsc(fn())
        decoy_strategies[(key, True, True)] = s.keepsc(s.keepsn(fn()))
        decoy_strategies[(key, 'M', True)] = s.keepsc(s.keepsm(fn()))

    for key, fn in _ENZYME_STRATEGIES:
        cleavage_strategies[(key, False, False)] = fn()
        cleavage_strategies[(key, True, False)] = s.keepsn(fn())
        cleavage_strategies[(key, 'M', False)] = s.keepsm(fn())

        cleavage_strategies[(key, False, True)] = s.keepsc(fn())
        cleavage_strategies[(key, True, True)] = s.keepsc(s.keepsn(fn()))
        cleavage_strategies[(key, 'M', True)] = s.keepsc(s.keepsm(fn()))


_register_all()

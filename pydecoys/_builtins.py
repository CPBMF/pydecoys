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

import itertools
from typing import Iterable, Literal, cast

from pydecoys import strategies as s


def reverse[T: s.SeqLike](sequence: T) -> T:
    """Return the reversed `sequence`."""
    return cast(T, sequence[::-1])


reverse_keepn = s.keepsn(reverse)
"""Return the reversed `sequence`, except N-terminal aa."""


reverse_keepc = s.keepsc(reverse)
"""Return the reversed `sequence`, except C-terminal aa."""


reverse_keepterm = s.keepsterm(reverse)
"""Return the reversed `sequence`, except terminal aas."""


def shuffle[T: s.SeqLike](sequence: T) -> T:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, "".join(new))


shuffle_keepn = s.keepsn(shuffle)
"""Return the shuffled `sequence`, except N-terminal aa."""


shuffle_keepc = s.keepsc(shuffle)
"""Return the shuffled `sequence`, except N-terminal aa."""


shuffle_keepterm = s.keepsterm(shuffle)
"""Return the shuffled `sequence`, except terminal aas."""


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


randomize = Randomize()
randomize_keepn = s.keepsn(Randomize())
randomize_keepc = s.keepsc(Randomize())
randomize_keepterm = s.keepsterm(Randomize())


decoy_strategy: dict[str, s.DecoyGenerator] = {
    "reverse": reverse,
    "reverse-keepn": reverse_keepn,
    "reverse-keepc": reverse_keepc,
    "reverse-keepterm": reverse_keepterm,
    "shuffle": shuffle,
    "shuffle-keepn": shuffle_keepn,
    "shuffle-keepc": shuffle_keepc,
    "shuffle-keepterm": shuffle_keepterm,
    "randomize": randomize,
    "randomize-keepn": randomize_keepn,
    "randomize-keepc": randomize_keepc,
    "randomize-keepterm": randomize_keepterm,
}


# Pre-defined pseudo-reverse and pseudo-shuffle DecoyGenerators

type _Enzyme = tuple[str, str | None, str | None, Literal['N', 'C'], str]

_TRYPSIN: _Enzyme = ('KR', 'P', None, 'C', 'trypsin')
_TRYPSIN_P: _Enzyme = ('KR', None, None, 'C', 'trypsinp')
_ARG_C: _Enzyme = ('R', 'P', None, 'C', 'argc')
_ASP_N: _Enzyme = ('BD', None, None, 'N', 'aspn')
_CHYMO: _Enzyme = ('FYWL', 'P', None, 'C', 'chymo')
_V8_DE: _Enzyme = ('BDEZ', 'P', None, 'C', 'v8de')
_LYS_C: _Enzyme = ('K', 'P', None, 'C', 'lysc')
_LYS_N: _Enzyme = ('K', None, None, 'N', 'lysn')
_PEPSIN_A: _Enzyme = ('FL', None, None, 'C', 'pepsina')
_CNBR: _Enzyme = ('M', None, None, 'C', 'cnbr')


def _register_enzymatic_strategies():
    strategies = itertools.product(
        [s.ReversePep, s.ShufflePep, s.RandomizePep],
        [
            _TRYPSIN,
            _TRYPSIN_P,
            _ARG_C,
            _ASP_N,
            _CHYMO,
            _V8_DE,
            _LYS_C,
            _LYS_N,
            _PEPSIN_A,
            _CNBR
        ],
        ['keepsn', 'keepsc', 'keepsterm', None],
    )

    for factory, enzyme, term in strategies:
        key = factory.__name__.lower() + '-' + enzyme[4]
        generator = factory(enzyme[0], enzyme[1], enzyme[2], enzyme[3])
        match term:
            case 'keepsn':
                generator = s.keepsn(generator)
                key += '-keepn'
            case 'keepsc':
                generator = s.keepsc(generator)
                key += '-keepc'
            case 'keepsterm':
                generator = s.keepsterm(generator)
                key += '-keepterm'
        decoy_strategy[key] = generator


_register_enzymatic_strategies()

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


# Pre-defined enzymes
# Those specifications were taken from <https://github.com/HUPO-PSI/psi-ms-CV>

type _Enzyme = tuple[str, str | None, str | None, Literal['N', 'C', 'both'], str]

_ALPHA_LP: _Enzyme = ('TASV', None, None, 'C', 'alphalp')
_ARG_C: _Enzyme = ('R', 'P', None, 'C', 'argc')
_ASP_N: _Enzyme = ('BD', None, None, 'N', 'aspn')
_ASP_N_AMBIC: _Enzyme = ('DE', None, None, 'N', 'aspnambic')
_CHYMO: _Enzyme = ('FYWL', 'P', None, 'C', 'chymo')
_CNBR: _Enzyme = ('M', None, None, 'C', 'cnbr')
_FORMIC_ACID: _Enzyme = ('D', None, None, 'both', 'formicacid')
_GLUTAMYL_ENDOPEPTIDASE: _Enzyme = ('E', None, 'E', 'C', 'gluc')
_LEUKOCYTE_ELASTASE: _Enzyme = ('ALIV', 'P', None, 'C', 'elastase')
_LYS_C: _Enzyme = ('K', 'P', None, 'C', 'lysc')
_LYS_C_P: _Enzyme = ('K', None, None, 'C', 'lyscp')
_LYS_N: _Enzyme = ('K', None, None, 'N', 'lysn')
_PEPSIN_A: _Enzyme = ('FL', None, None, 'C', 'pepsina')
_TRYPSIN: _Enzyme = ('KR', 'P', None, 'C', 'trypsin')
_TRYPSIN_P: _Enzyme = ('KR', None, None, 'C', 'trypsinp')
_TRYP_CHYMO: _Enzyme = ('FYWLKR', 'P', None, 'C', 'trypchymo')
_TRYP_N: _Enzyme = ('KR', None, None, 'N', 'trypn')
_TWO_IODOBENZOATE: _Enzyme = ('W', None, None, 'C', '2iodobenzoate')
_V8_DE: _Enzyme = ('BDEZ', 'P', None, 'C', 'v8de')
_V8_E: _Enzyme = ('EZ', 'P', None, 'C', 'v8e')

# Pro-C is the only one that needs a regex pattern
_PROLINE_ENDOPEPTIDASE = (r'([HKR]P)(?!P)', 'proc')


def _register_enzymatic_strategies():
    strategies = itertools.product(
        [s.ReversePep, s.ShufflePep, s.RandomizePep],
        [
            _ALPHA_LP,
            _ARG_C,
            _ASP_N,
            _ASP_N_AMBIC,
            _CHYMO,
            _CNBR,
            _FORMIC_ACID,
            _GLUTAMYL_ENDOPEPTIDASE,
            _LEUKOCYTE_ELASTASE,
            _LYS_C,
            _LYS_C_P,
            _LYS_N,
            _PEPSIN_A,
            _TRYPSIN,
            _TRYPSIN_P,
            _TRYP_CHYMO,
            _TRYP_N,
            _TWO_IODOBENZOATE,
            _V8_DE,
            _V8_E,
            None,
        ],
        ['keepsn', 'keepsc', 'keepsterm', None],
    )

    # TODO: Find a better way of coding this
    # Maybe exposing a function that automatically registers all combinations
    # of an enzyme/EnzymeSpecificGenerator might be useful for plug-ins?
    # Not make it automatic so plug-ins can opt out of it?
    for factory, enzyme, term in strategies:
        if enzyme is None:
            key = factory.__name__.lower() + '-' + _PROLINE_ENDOPEPTIDASE[1]
            generator = factory.from_regex(_PROLINE_ENDOPEPTIDASE[0])
        else:
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

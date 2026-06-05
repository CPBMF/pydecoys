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

from typing import Callable, Iterable, Literal, cast

from pydecoys import strategies as s


decoy_strategy: dict[str, s.DecoyGenerator] = {}


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


_PLAIN_STRATEGIES: list[tuple[str, Callable[[], s.DecoyGenerator]]] = [
    ('reverse',   lambda: reverse),  # noqa: E272
    ('shuffle',   lambda: shuffle),  # noqa: E272
    ('randomize', Randomize)
]


_ENZYME_STRATEGIES: list[tuple[str, type[s.EnzymeSpecificGenerator]]] = [
    ('reversepep',   s.ReversePep),
    ('shufflepep',   s.ShufflePep),
    ('randomizepep', s.RandomizePep)
]


_TERMINAL_FACTORIES: list[
    tuple[str, Callable[[s.DecoyGenerator], s.DecoyGenerator]]
] = [
    ('keepn',    s.keepsn),
    ('keepc',    s.keepsc),
    ('keepterm', s.keepsterm)
]


def _register_all() -> None:
    for base_key, make in _PLAIN_STRATEGIES:
        decoy_strategy[base_key] = make()
        for suffix, factory in _TERMINAL_FACTORIES:
            decoy_strategy[f'{base_key}-{suffix}'] = factory(make())

    for class_key, cls in _ENZYME_STRATEGIES:
        for regex, sense, enzyme_key in _ENZYMES:
            base_key = f'{class_key}-{enzyme_key}'
            decoy_strategy[base_key] = cls(regex, sense)

            for suffix, factory in _TERMINAL_FACTORIES:
                decoy_strategy[f'{base_key}-{suffix}'] = factory(cls(regex, sense))


_register_all()

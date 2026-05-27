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
    _AA_TO_INDEX = {aa: i for i, aa in enumerate(s.AMINOACIDS)}

    def __init__(self):
        self._weights = [0] * 20

    def learn_context(self, sequences: Iterable[s.SeqLike]):
        self._weights = [0] * 20

        for seq in sequences:
            for aa in seq:
                idx = self._AA_TO_INDEX.get(aa)
                if idx is not None:
                    self._weights[idx] += 1

    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        length = len(sequence)
        new = s.RAND.choices(s.AMINOACIDS, weights=self._weights, k=length)
        return s.seq_cast(sequence, "".join(new))


randomize = Randomize()
randomize_keepn = s.keepsn(Randomize())
randomize_keepc = s.keepsc(Randomize())
randomize_keepterm = s.keepsterm(Randomize())


# Pre-defined pseudo-reverse and pseudo-shuffle DecoyGenerators

type _Enzyme = tuple[str, str | None, Literal['N', 'C']]

_TRYPSIN: _Enzyme = ('KR', 'P', 'C')
_STRICT_TRYPSIN: _Enzyme = ('KR', None, 'C')
_ARG_C: _Enzyme = ('R', 'P', 'C')
_ASP_N: _Enzyme = ('D', None, 'N')
_CHYMO: _Enzyme = ('FLWY', 'P', 'C')
_GLU_C: _Enzyme = ('DE', 'P', 'C')
_LYS_C: _Enzyme = ('K', 'P', 'C')
_LYS_N: _Enzyme = ('K', None, 'N')
_PEPSIN_A: _Enzyme = ('FL', None, 'C')
_CNBR: _Enzyme = ('M', None, 'C')


def _with_enzyme[T: s.EnzymeSpecificGenerator](
    generator: type[T],
    enzyme: _Enzyme,
) -> T:
    return generator(enzyme[0], enzyme[1], enzyme[2])


reversepep_trypsin = _with_enzyme(s.ReversePep, _TRYPSIN)
reversepep_stricttrypsin = _with_enzyme(s.ReversePep, _STRICT_TRYPSIN)
reversepep_argc = _with_enzyme(s.ReversePep, _ARG_C)
reversepep_aspn = _with_enzyme(s.ReversePep, _ASP_N)
reversepep_chymo = _with_enzyme(s.ReversePep, _CHYMO)
reversepep_gluc = _with_enzyme(s.ReversePep, _GLU_C)
reversepep_lysc = _with_enzyme(s.ReversePep, _LYS_C)
reversepep_lysn = _with_enzyme(s.ReversePep, _LYS_N)
reversepep_pepsina = _with_enzyme(s.ReversePep, _PEPSIN_A)
reversepep_cnbr = _with_enzyme(s.ReversePep, _CNBR)

# Keep-n versions
reversepep_stricttrypsin_keepn = \
    s.keepsn(_with_enzyme(s.ReversePep, _STRICT_TRYPSIN))

shufflepep_trypsin = _with_enzyme(s.ShufflePep, _TRYPSIN)
shufflepep_stricttrypsin = _with_enzyme(s.ShufflePep, _STRICT_TRYPSIN)
shufflepep_argc = _with_enzyme(s.ShufflePep, _ARG_C)
shufflepep_aspn = _with_enzyme(s.ShufflePep, _ASP_N)
shufflepep_chymo = _with_enzyme(s.ShufflePep, _CHYMO)
shufflepep_gluc = _with_enzyme(s.ShufflePep, _GLU_C)
shufflepep_lysc = _with_enzyme(s.ShufflePep, _LYS_C)
shufflepep_lysn = _with_enzyme(s.ShufflePep, _LYS_N)
shufflepep_pepsina = _with_enzyme(s.ShufflePep, _PEPSIN_A)
shufflepep_cnbr = _with_enzyme(s.ShufflePep, _CNBR)

# Keep-n versions
shufflepep_stricttrypsin_keepn = \
    s.keepsn(_with_enzyme(s.ShufflePep, _STRICT_TRYPSIN))


randomizepep_trypsin = _with_enzyme(s.RandomizePep, _TRYPSIN)
randomizepep_stricttrypsin = _with_enzyme(s.RandomizePep, _STRICT_TRYPSIN)
randomizepep_argc = _with_enzyme(s.RandomizePep, _ARG_C)
randomizepep_aspn = _with_enzyme(s.RandomizePep, _ASP_N)
randomizepep_chymo = _with_enzyme(s.RandomizePep, _CHYMO)
randomizepep_gluc = _with_enzyme(s.RandomizePep, _GLU_C)
randomizepep_lysc = _with_enzyme(s.RandomizePep, _LYS_C)
randomizepep_lysn = _with_enzyme(s.RandomizePep, _LYS_N)
randomizepep_pepsina = _with_enzyme(s.RandomizePep, _PEPSIN_A)
randomizepep_cnbr = _with_enzyme(s.RandomizePep, _CNBR)

randomizepep_stricttrypsin_keepn = \
    s.keepsn(_with_enzyme(s.RandomizePep, _STRICT_TRYPSIN))


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
    "reversepep-trypsin": reversepep_trypsin,
    "reversepep-stricttrypsin": reversepep_stricttrypsin,
    "reversepep-argc": reversepep_argc,
    "reversepep-aspn": reversepep_aspn,
    "reversepep-chymo": reversepep_chymo,
    "reversepep-gluc": reversepep_gluc,
    "reversepep-lysc": reversepep_lysc,
    "reversepep-lysn": reversepep_lysn,
    "reversepep-stricttrypsin-keepn": reversepep_stricttrypsin_keepn,  # noqa: E501
    "shufflepep-trypsin": shufflepep_trypsin,
    "shufflepep-stricttrypsin": shufflepep_stricttrypsin,
    "shufflepep-argc": shufflepep_argc,
    "shufflepep-aspn": shufflepep_aspn,
    "shufflepep-chymo": shufflepep_chymo,
    "shufflepep-gluc": shufflepep_gluc,
    "shufflepep-lysc": shufflepep_lysc,
    "shufflepep-lysn": shufflepep_lysn,
    "shufflepep-stricttrypsin-keepn": shufflepep_stricttrypsin_keepn,  # noqa: E501
    "randomizepep-trypsin": randomizepep_trypsin,
    "randomizepep-stricttrypsin": randomizepep_stricttrypsin,
    "randomizepep-argc": randomizepep_argc,
    "randomizepep-aspn": randomizepep_aspn,
    "randomizepep-chymo": randomizepep_chymo,
    "randomizepep-gluc": randomizepep_gluc,
    "randomizepep-lysc": randomizepep_lysc,
    "randomizepep-lysn": randomizepep_lysn,
    "randomizepep-stricttrypsin-keepn": randomizepep_stricttrypsin_keepn,  # noqa: E501
}

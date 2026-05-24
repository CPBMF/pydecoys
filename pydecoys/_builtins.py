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

from __future__ import annotations

from typing import Literal, overload, TYPE_CHECKING

if TYPE_CHECKING:
    from Bio.Seq import Seq, MutableSeq

from pydecoys import strategies as s

# This file has a lot of overloads, but sphinx autodoc doesn't seem to catch
# the overloads if I move them to a .pyi file.


@overload
def reverse(sequence: Seq) -> Seq: ...


@overload
def reverse(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse(sequence: str) -> str: ...


def reverse(sequence: s.SeqLike) -> s.SeqLike:
    """Return the reversed `sequence`."""
    return sequence[::-1]


@overload
def reverse_keep_n(sequence: Seq) -> Seq: ...


@overload
def reverse_keep_n(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse_keep_n(sequence: str) -> str: ...


def reverse_keep_n(sequence: s.SeqLike) -> s.SeqLike:
    """Return the reversed `sequence`, except N-terminal aa."""
    return sequence[0] + sequence[:0:-1]


@overload
def reverse_keep_c(sequence: Seq) -> Seq: ...


@overload
def reverse_keep_c(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse_keep_c(sequence: str) -> str: ...


def reverse_keep_c(sequence: s.SeqLike) -> s.SeqLike:
    """Return the reversed `sequence`, except C-terminal aa."""
    return sequence[-2::-1] + sequence[-1]


@overload
def reverse_keep_term(sequence: Seq) -> Seq: ...


@overload
def reverse_keep_term(sequence: MutableSeq) -> MutableSeq: ...


@overload
def reverse_keep_term(sequence: str) -> str: ...


def reverse_keep_term(sequence: s.SeqLike) -> s.SeqLike:
    """Return the reversed `sequence`, except terminal aas."""
    return sequence[0] + sequence[-2:0:-1] + sequence[-1]


@overload
def shuffle(sequence: Seq) -> Seq: ...


@overload
def shuffle(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle(sequence: str) -> str: ...


def shuffle(sequence: s.SeqLike) -> s.SeqLike:
    """Return the shuffled `sequence`."""
    new = list(sequence)
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, "".join(new))


@overload
def shuffle_keep_n(sequence: Seq) -> Seq: ...


@overload
def shuffle_keep_n(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle_keep_n(sequence: str) -> str: ...


def shuffle_keep_n(sequence: s.SeqLike) -> s.SeqLike:
    """Return the shuffled `sequence`, except N-terminal aa."""
    new = list(sequence[1:])
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, sequence[0] + "".join(new))


@overload
def shuffle_keep_c(sequence: Seq) -> Seq: ...


@overload
def shuffle_keep_c(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle_keep_c(sequence: str) -> str: ...


def shuffle_keep_c(sequence: s.SeqLike) -> s.SeqLike:
    """Return the shuffled `sequence`, except C-terminal aa."""
    new = list(sequence[:-1])
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, "".join(new) + sequence[-1])


@overload
def shuffle_keep_term(sequence: Seq) -> Seq: ...


@overload
def shuffle_keep_term(sequence: MutableSeq) -> MutableSeq: ...


@overload
def shuffle_keep_term(sequence: str) -> str: ...


def shuffle_keep_term(sequence: s.SeqLike) -> s.SeqLike:
    """Return the shuffled `sequence`, except terminal aas."""
    new = list(sequence[1:-1])
    s.RAND.shuffle(new)
    return s.seq_cast(sequence, sequence[0] + "".join(new) + sequence[-1])


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
    keep_term: Literal['N', 'C', 'both', None] = None,
) -> T:
    return generator(enzyme[0], enzyme[1], enzyme[2], keep_term)


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
reversepep_stricttrypsin_keepn = _with_enzyme(s.ReversePep, _STRICT_TRYPSIN, 'N')

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
shufflepep_stricttrypsin_keepn = _with_enzyme(s.ShufflePep, _STRICT_TRYPSIN, 'N')

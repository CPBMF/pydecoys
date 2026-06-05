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

"""Decoy-generation logic, as well as API to write new decoy strategies.

The main API of `strategies` is the :type:`DecoyGenerator` type. This is a
simple type alias that takes a type ``T: SeqLike`` target and returns a type
``T`` decoy. The :type:`SeqLike` is an alias for ``'str | Seq | MutableSeq'``.

For decoy strategies that need context from the target database (for example,
that use a Markov State Model), implement the :class:`ContextfulGenerator`
protocol.

The :class:`ReversePep`, :class:`ShufflePep` and :class:`RandomizePep` classes
allow easy definition of new enzyme specifications for pseudo-decoy
strategies via instantiation. To define new strategies using enzymatic
peptides, subclass from :class:`EnzymeSpecificGenerator`.
"""

__all__ = [
    'RAND',
    'STD_AMINOACIDS',
    'EXT_AMINOACIDS',
    'SeqLike',
    'DecoyGenerator',
    'ContextfulGenerator',
    'EnzymeSpecificGenerator',
    'RandomizePep',
    'ReversePep',
    'ShufflePep',
    'seq_cast',
    'keepsn',
    'keepsc',
    'keepsterm',
    'fuses',
]


from pydecoys.strategies.core import (
    RAND,
    STD_AMINOACIDS,
    EXT_AMINOACIDS,
    SeqLike,
    DecoyGenerator,
    ContextfulGenerator,
    EnzymeSpecificGenerator,
    RandomizePep,
    ReversePep,
    ShufflePep,
    seq_cast
)
from pydecoys.strategies.factories import keepsc, keepsn, keepsterm, fuses

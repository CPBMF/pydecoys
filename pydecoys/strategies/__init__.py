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

The main API of `strategies` is the :class:`DecoyGenerator` type. This is a
simple protocol that only implements a `__call__` function and appropriate
type overloads.

For decoy strategies that need context from the target database (for example,
that use a Markov State Model), implement the :class:`ContextfulGenerator`
protocol.

The :class:`ReversePep` and :class:`ShufflePep` classes allow easy definition
of new enzyme specifications for pseudo-reverse and pseudo-shuffle strategies
via instantiation. To define new strategies using enzymatic peptides, subclass
from :class:`EnzymeSpecificGenerator`.
"""


from .core import (
    SeqLike,
    Seq_,
    MutableSeq_,
    DecoyGenerator,
    ContextfulGenerator,
    EnzymeSpecificGenerator,
    ReversePep,
    ShufflePep,
    RAND,
    seq_cast,
    AMINOACIDS,
)


__all__ = [
    'SeqLike',
    'Seq_',
    'MutableSeq_',
    'DecoyGenerator',
    'ContextfulGenerator',
    'EnzymeSpecificGenerator',
    'ReversePep',
    'ShufflePep',
    'RAND',
    'seq_cast',
    'AMINOACIDS',
]

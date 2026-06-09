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

"""Public API to define and register new decoy strategies.

The main API of `strategies` is the :type:`DecoyGenerator` type. This is a
simple type alias that takes a type ``T: SeqLike`` target and returns a type
``T`` decoy. The :type:`SeqLike` is an alias for ``'str | Seq | MutableSeq'``.

For decoy strategies that need context from the target database (for example,
that use a Markov State Model), implement the :class:`ContextfulGenerator`
protocol.

To define new strategies using enzymatic peptides, subclass from
:class:`EnzymeSpecificGenerator`.

To register strategies, three decorators are available:
:func:`register_function`, :func:`register_class` and
:func:`register_cleavage_aware`.

Available functionality
-----------------------
:data:`RAND`
    Pseudo-Random Number Generator with a fixed seed for stochastic strategies.
:data:`STD_AMINOACIDS`
    String of the standard 20 proteinogenic amino acids single-letter codes.
:data:`EXT_AMINOACIDS`
    String of the 22 proteinogenic amino acids and three special single-letter
    codes.
:type:`DecoyGenerator`
    Type alias defining the type of a decoy strategy function.
:class:`ContextfulGenerator`
    Protocol for decoy strategies that require previous knowledge of the
    target database.
:class:`EnzymeSpecificGenerator`
    ABC for decoy strategies create decoy cleavaged peptides.
:func:`keepsn`
    Return a modified strategy that doesn't alter the N-terminal aa.
:func:`keepsc`
    Return a modified strategy that doesn't alter the C-terminal aa.
:func:`keepsterm`
    Return a modified strategy that doesn't alter the terminal aas.
:func:`register_function`
    Register a string key for the decorated function.
:func:`register_class`
    Register a string key for the decorated callable class.
:func:`register_cleavage_aware`
    Register string keys with each cleavage agent for the decorated class.
:func:`add_cleavage_agent`
    Add a new cleavage agent and string key to cleavage agent registry.
:func:`add_callable`
    Add a new callable and string key to strategies registry.
:func:`view_strategies`
    Return a read-only view of all strategy string keys.
:func:`view_cleavage_agents`
    Return a read-only view of the dictionary of cleavage agents.
"""

__all__ = [
    'RAND',
    'STD_AMINOACIDS',
    'EXT_AMINOACIDS',
    'SeqLike',
    'DecoyGenerator',
    'ContextfulGenerator',
    'EnzymeSpecificGenerator',
    'seq_cast',
    'keepsn',
    'keepsc',
    'keepsterm',
    'register_class',
    'register_cleavage_aware',
    'register_function',
    'register_callable',
    'register_cleavage_agent',
    'view_cleavage_agents',
    'view_strategies',
]


from pydecoys.strategies.core import (
    RAND,
    STD_AMINOACIDS,
    EXT_AMINOACIDS,
    SeqLike,
    DecoyGenerator,
    ContextfulGenerator,
    EnzymeSpecificGenerator,
    seq_cast
)
from pydecoys.strategies.factories import keepsc, keepsn, keepsterm
from pydecoys.strategies.registry import (
    register_class,
    register_cleavage_aware,
    register_function,
    register_callable,
    register_cleavage_agent,
    view_cleavage_agents,
    view_strategies,
)

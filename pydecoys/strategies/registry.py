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

"""Registry for strategies and related functionality."""

from collections.abc import KeysView
from types import MappingProxyType
from typing import Any, Callable, Literal

from pydecoys import strategies as s


type _Enzyme = tuple[str, Literal['N', 'C', 'both']]


DECOY_STRATEGIES: dict[str, s.DecoyGenerator] = {}
CLEAVAGE_AGENTS: dict[str, _Enzyme] = {}
CLEAVAGE_STRATEGIES: list[
    tuple[
        str,
        type[s.EnzymeSpecificGenerator],
        bool,
        bool,
        bool,
        dict[str, Any]
    ]
] = []


def register_function(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True,
):
    """Register the decorated function under the given str key.

    The str key will become accessible as an available strategy in
    :mod:`pydecoys`, mapping to the function.

    Parameters
    ----------
    str_key
        Lower case string identifying the decoy strategy. Must not be already
        defined.
    keep_n
        If `True`, a version of the function that doesn't alter the N-terminal
        aminoacid will be registered automatically. The key will be
        ``f'{str_key}-keepn'``.
    keep_c
        If `True`, a version of the function that doesn't alter the C-terminal
        aminoacid will be registered automatically. The key will be
        ``f'{str_key}-keepc'``.
    keep_term
        If `True`, a version of the function that doesn't alter the terminal
        aminoacids will be registered automatically. The key will be
        ``f'{str_key}-keepterm'``.
    """
    def decorator[T: s.DecoyGenerator](fn: T) -> T:
        if keep_n:
            register_callable(f'{str_key}-keepn', s.keepsn(fn))
        if keep_c:
            register_callable(f'{str_key}-keepc', s.keepsc(fn))
        if keep_term:
            register_callable(f'{str_key}-keepterm', s.keepsterm(fn))

        register_callable(str_key, fn)

        return fn
    return decorator


def register_class(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True,
    **kwargs: Any,
):
    """Register the decorated class under the given str key.

    The str key will become accessible as an available strategy in
    :mod:`pydecoys`, mapping to an instance of the class.

    Parameters
    ----------
    str_key
        Lower case string identifying the decoy strategy. Must not be already
        defined.
    keep_n
        If `True`, an instance of the class that doesn't alter the N-terminal
        aminoacid will be registered automatically. The key will be
        ``f'{str_key}-keepn'``.
    keep_c
        If `True`, an instance of the class that doesn't alter the C-terminal
        aminoacid will be registered automatically. The key will be
        ``f'{str_key}-keepc'``.
    keep_term
        If `True`, an instance of the class that doesn't alter the terminal
        aminoacids will be registered automatically. The key will be
        ``f'{str_key}-keepterm'``.
    kwargs
        Optional keyword arguments that'll be passed to the class during
        initialization.
    """
    def decorator[C: Callable[[], s.DecoyGenerator]](cls: C) -> C:
        if keep_n:
            register_callable(f'{str_key}-keepn', s.keepsn(cls(**kwargs)))
        if keep_c:
            register_callable(f'{str_key}-keepc', s.keepsc(cls(**kwargs)))
        if keep_term:
            register_callable(f'{str_key}-keepterm', s.keepsterm(cls(**kwargs)))

        register_callable(str_key, cls(**kwargs))

        return cls
    return decorator


def register_cleavage_aware(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True,
    **kwargs: Any
):
    """Register the decorated class under the given str key for each available
    cleavage agent.

    The combinations of ``f'{str_key}-<enzyme>'`` will become accessible as
    available strategies in :mod:`pydecoys`, mapping to an instance of the
    class. If a new cleavage agent is added later, a new instance of the class
    will be created and registered for it.

    Parameters
    ----------
    str_key
        Lower case string identifying the decoy strategy. Must not be already
        defined.
    keep_n
        If `True`, an instance of the class that doesn't alter the N-terminal
        aminoacid will be registered automatically. The key will be
        ``f'{str_key}-keepn'``.
    keep_c
        If `True`, an instance of the class that doesn't alter the C-terminal
        aminoacid will be registered automatically. The key will be
        ``f'{str_key}-keepc'``.
    keep_term
        If `True`, an instance of the class that doesn't alter the terminal
        aminoacids will be registered automatically. The key will be
        ``f'{str_key}-keepterm'``.
    kwargs
        Optional keyword arguments that'll be passed to the class during
        initialization.
    """
    def decorator(cls: type[s.EnzymeSpecificGenerator]):
        CLEAVAGE_STRATEGIES.append((
            str_key,
            cls,
            keep_n,
            keep_c,
            keep_term,
            kwargs
        ))

        for agent, (regex, sense) in CLEAVAGE_AGENTS.items():
            key = f'{str_key}-{agent}'
            if keep_n:
                register_callable(
                    f'{key}-keepn',
                    s.keepsn(cls(regex, sense, **kwargs))
                )
            if keep_c:
                register_callable(
                    f'{key}-keepc',
                    s.keepsc(cls(regex, sense, **kwargs))
                )
            if keep_term:
                register_callable(
                    f'{key}-keepterm',
                    s.keepsterm(cls(regex, sense, **kwargs))
                )

            register_callable(key, cls(regex, sense, **kwargs))
        return cls
    return decorator


def register_cleavage_agent(
    str_key: str,
    pattern: str,
    sense: Literal['N', 'C', 'both']
):
    """Add the specified cleavage agent under the given str key.

    All current and added later :class:`EnzymeSpecificGenerator` strategies
    will have an instance registered as ``f'<strategy>-{str_key}'``.

    Parameters
    ----------
    str_key
        Lower case string identifying the cleavage agent. Must not be already
        defined.
    pattern
        String for the cleavage agent's regex pattern. The regex MUST match
        **only** the cleavage sites that shouldn't be altered. The cleavage
        sites MUST be captured by the regex pattern.
    sense
        Whether the enzyme cleaves the C-terminal, N-terminal or both termini
        of the cleavage site.

    Examples
    --------
    >>> add_cleavage_agent('agent', r'(K)', 'C')
    >>> seq = 'DNIDYKAVYR'
    >>> seq_as_decoy(seq, 'reversepep-agent')
    'YDINDKRYVA'

    The key must be a lowercase string:

    >>> add_cleavage_agent(5, random_seq)
    Traceback (most recent call last):
        ...
    TypeError: Need a string for the key (lower case)
    >>> add_cleavage_agent('AGENT', r'(K)', 'C')
    Traceback (most recent call last):
        ...
    ValueError: Key 'AGENT' should be lower case

    It must not be already defined:

    >>> add_callable('agent', r'(K)', 'C')
    Traceback (most recent call last):
        ...
    ValueError: Key 'agent' already defined
    """

    _validate_key(str_key)

    if str_key in CLEAVAGE_AGENTS:
        raise ValueError(f"Key '{str_key}' already defined")

    CLEAVAGE_AGENTS[str_key] = (pattern, sense)

    for base_key, fn, keep_n, keep_c, keep_term, kwargs in CLEAVAGE_STRATEGIES:
        key = f'{base_key}-{str_key}'
        if keep_n:
            register_callable(f'{key}-keepn', s.keepsn(fn(pattern, sense, **kwargs)))
        if keep_c:
            register_callable(f'{key}-keepc', s.keepsc(fn(pattern, sense, **kwargs)))
        if keep_term:
            register_callable(
                f'{key}-keepterm',
                s.keepsterm(fn(pattern, sense, **kwargs))
            )
        register_callable(f'{key}', fn(pattern, sense, **kwargs))


def register_callable(strategy_key: str, strategy_fn: s.DecoyGenerator):
    """Add a new decoy strategy to the strategy registry.

    The callable will be registered as-is, without other side-effects.

    Parameters
    ----------
    strategy_key
        Lower case string identifying the decoy strategy. Must not be already
        defined.
    strategy_fn
        A function following the :type:`strategies.DecoyGenerator` signature.

    Examples
    --------
    Given a `random_seq` function that takes a sequence and returns a new,
    unrelated sequence of same size:

    >>> def random_seq(sequence): ...
    >>> add_callable('randomseq', random_seq)
    >>> seq = 'DNIDYKAVYR'
    >>> seq_as_decoy(seq, 'randomseq')  # doctest: +SKIP
    'LLEETLSWQC'

    The strategy key must be a lowercase string:

    >>> add_callable(5, random_seq)
    Traceback (most recent call last):
        ...
    TypeError: Need a string for the key (lower case)
    >>> add_callable('RANDOMSEQ', random_seq)
    Traceback (most recent call last):
        ...
    ValueError: Key 'RANDOMSEQ' should be lower case

    It must not be already defined:

    >>> add_callable('randomseq', random_seq)
    Traceback (most recent call last):
        ...
    ValueError: Key 'randomseq' already defined

    The second argument must be a callable:

    >>> add_callable('randomseq2', 'randomseq2')
    Traceback (most recent call last):
        ...
    TypeError: Strategy function must be a callable
    """

    _validate_key(strategy_key)

    if strategy_key in DECOY_STRATEGIES:
        raise ValueError(f"Key '{strategy_key}' already defined")
    if not callable(strategy_fn):
        raise TypeError("Strategy function must be a callable")

    DECOY_STRATEGIES[strategy_key] = strategy_fn


def view_strategies() -> KeysView[str]:
    return DECOY_STRATEGIES.keys()


def view_cleavage_agents() -> MappingProxyType[str, _Enzyme]:
    return MappingProxyType(CLEAVAGE_AGENTS)


def _validate_key(key):
    if not isinstance(key, str):
        raise TypeError("Need a string for the key (lower case)")
    if not key:
        raise ValueError("Key required (lower case string)")
    if not key.islower():
        raise ValueError(f"Key '{key}' should be lower case")

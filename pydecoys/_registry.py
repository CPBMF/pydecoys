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


def register_fn(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True,
):
    def decorator[T: s.DecoyGenerator](fn: T) -> T:
        if keep_n:
            add_callable(f'{str_key}-keepn', s.keepsn(fn))
        if keep_c:
            add_callable(f'{str_key}-keepc', s.keepsc(fn))
        if keep_term:
            add_callable(f'{str_key}-keepterm', s.keepsterm(fn))

        add_callable(str_key, fn)

        return fn
    return decorator


def register_cls(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True,
    **kwargs: Any,
):
    def decorator[C: Callable[[], s.DecoyGenerator]](cls: C) -> C:
        if keep_n:
            add_callable(f'{str_key}-keepn', s.keepsn(cls(**kwargs)))
        if keep_c:
            add_callable(f'{str_key}-keepc', s.keepsc(cls(**kwargs)))
        if keep_term:
            add_callable(f'{str_key}-keepterm', s.keepsterm(cls(**kwargs)))

        add_callable(str_key, cls(**kwargs))

        return cls
    return decorator


def register_cleavage_aware(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True,
    **kwargs: Any
):
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
                add_callable(
                    f'{key}-keepn',
                    s.keepsn(cls(regex, sense, **kwargs))
                )
            if keep_c:
                add_callable(
                    f'{key}-keepc',
                    s.keepsc(cls(regex, sense, **kwargs))
                )
            if keep_term:
                add_callable(
                    f'{key}-keepterm',
                    s.keepsterm(cls(regex, sense, **kwargs))
                )

            add_callable(key, cls(regex, sense, **kwargs))
        return cls
    return decorator


def add_cleavage_agent(str_key: str, pattern: str, sense: Literal['N', 'C', 'both']):
    CLEAVAGE_AGENTS[str_key] = (pattern, sense)
    for base_key, fn, keep_n, keep_c, keep_term, kwargs in CLEAVAGE_STRATEGIES:
        key = f'{base_key}-{str_key}'
        if keep_n:
            add_callable(f'{key}-keepn', s.keepsn(fn(pattern, sense, **kwargs)))
        if keep_c:
            add_callable(f'{key}-keepc', s.keepsc(fn(pattern, sense, **kwargs)))
        if keep_term:
            add_callable(f'{key}-keepterm', s.keepsterm(fn(pattern, sense, **kwargs)))
        add_callable(f'{key}', fn(pattern, sense, **kwargs))


def add_callable(strategy_key: str, strategy_fn: s.DecoyGenerator):
    """Add a new decoy strategy.

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
    TypeError: Need a string for the decoy strategy (lower case)
    >>> add_callable('RANDOMSEQ', random_seq)
    Traceback (most recent call last):
        ...
    ValueError: Strategy key 'RANDOMSEQ' should be lower case

    It must not be already defined:

    >>> add_callable('randomseq', random_seq)
    Traceback (most recent call last):
        ...
    ValueError: Strategy key 'randomseq' already defined

    The second argument must be a callable:

    >>> add_callable('randomseq2', 'randomseq2')
    Traceback (most recent call last):
        ...
    TypeError: Strategy function must be a callable
    """

    if not isinstance(strategy_key, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy_key:
        raise ValueError("Strategy required (lower case string)")
    if not strategy_key.islower():
        raise ValueError(f"Strategy key '{strategy_key}' should be lower case")

    if strategy_key in DECOY_STRATEGIES:
        raise ValueError(f"Strategy key '{strategy_key}' already defined")
    DECOY_STRATEGIES[strategy_key] = strategy_fn


def view_strategies() -> KeysView[str]:
    return DECOY_STRATEGIES.keys()


def view_cleavage_agents() -> MappingProxyType[str, _Enzyme]:
    return MappingProxyType(CLEAVAGE_AGENTS)

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
    args: list[list[Any]] = []
):
    def decorator[C: Callable[[], s.DecoyGenerator]](cls: C) -> C:
        if keep_n:
            add_callable(f'{str_key}-keepn', s.keepsn(cls()))
        if keep_c:
            add_callable(f'{str_key}-keepc', s.keepsc(cls()))
        if keep_term:
            add_callable(f'{str_key}-keepterm', s.keepsterm(cls()))

        add_callable(str_key, cls())

        return cls
    return decorator


def register_cleavage_aware(
    str_key: str,
    keep_n=True,
    keep_c=True,
    keep_term=True
):
    def decorator(cls: type[s.EnzymeSpecificGenerator]):
        CLEAVAGE_STRATEGIES.append((str_key, cls, keep_n, keep_c, keep_term))

        for agent, (regex, sense) in CLEAVAGE_AGENTS.items():
            key = f'{str_key}-{agent}'
            if keep_n:
                add_callable(f'{key}-keepn', s.keepsn(cls(regex, sense)))
            if keep_c:
                add_callable(f'{key}-keepc', s.keepsc(cls(regex, sense)))
            if keep_term:
                add_callable(f'{key}-keepterm', s.keepsterm(cls(regex, sense)))

            add_callable(key, cls(regex, sense))
        return cls
    return decorator


def add_cleavage_agent(str_key: str, pattern: str, sense: Literal['N', 'C', 'both']):
    CLEAVAGE_AGENTS[str_key] = (pattern, sense)
    for base_key, fn, keep_n, keep_c, keep_term in CLEAVAGE_STRATEGIES:
        key = f'{base_key}-{str_key}'
        if keep_n:
            add_callable(f'{key}-keepn', s.keepsn(fn(pattern, sense)))
        if keep_c:
            add_callable(f'{key}-keepc', s.keepsc(fn(pattern, sense)))
        if keep_term:
            add_callable(f'{key}-keepterm', s.keepsterm(fn(pattern, sense)))
        add_callable(f'{key}', fn(pattern, sense))


def add_callable(str_key: str, fn):
    DECOY_STRATEGIES[str_key] = fn


def strategies() -> KeysView[str]:
    return DECOY_STRATEGIES.keys()


def cleavage_agents() -> MappingProxyType[str, _Enzyme]:
    return MappingProxyType(CLEAVAGE_AGENTS)

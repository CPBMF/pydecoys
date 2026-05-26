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

"""Factory functions to generate versions of DecoyGenerator that don't alter
terminal aminoacids.
"""

from functools import wraps
from typing import Iterable

from pydecoys.strategies import DecoyGenerator, SeqLike
from pydecoys.strategies.core import ContextfulGenerator


def keepsn[T: SeqLike](fn: DecoyGenerator[T]) -> DecoyGenerator[T]:
    """Factory that transform a DecoyGenerator into a new DecoyGenerator that
    doesn't alter the N-terminal aa.

    The `sequence` value is passed directly without the aminoacid that should
    be preserved. The aa is reinserted after the decorated function returns.

    :py:class:`strategies.ContextfulGenerator` objects preserve the
    :py:meth:`strategies.ContextfulGenerator.learn_context` method, but will
    receive each sequence without its N-terminal aa.

    Parameters
    ----------
    fn
        A :py:type:`strategies.DecoyGenerator`.

    Returns
    -------
    A version of `fn` that doesn't alter the N-terminal aminoacid of the target
    protein. If `fn` is a :py:class:`strategies.ContextfulGenerator`, the
    returned function will also be.

    Notes
    -----
    This function returns a new closure, meaning attributes other than
    metadata are lost.
    """
    @wraps(fn)
    def wrapper(sequence):
        return sequence[0] + fn(sequence[1:])

    if isinstance(fn, ContextfulGenerator):

        @wraps(fn.learn_context)
        def learn_context(sequences: Iterable[SeqLike]):
            sequences = (sequence[1:] for sequence in sequences)
            return fn.learn_context(sequences)
        wrapper.learn_context = learn_context  # type: ignore

    return wrapper


def keepsc[T: SeqLike](fn: DecoyGenerator[T]) -> DecoyGenerator[T]:
    """Decorator that transform a DecoyGenerator into a new DecoyGenerator that
    doesn't alter the C-terminal aa.

    The `sequence` value is passed directly without the aminoacid that should
    be preserved. The aa is reinserted after the decorated function returns.

    :py:class:`strategies.ContextfulGenerator` objects preserve the
    :py:meth:`strategies.ContextfulGenerator.learn_context` method, but will
    receive each sequence without its C-terminal aa.

    Parameters
    ----------
    fn
        A :py:type:`strategies.DecoyGenerator`.

    Returns
    -------
    A version of `fn` that doesn't alter the C-termina aminoacid of the target
    protein. If `fn` is a :py:class:`strategies.ContextfulGenerator`, the
    returned function will also be.

    Notes
    -----
    This function returns a new closure, meaning attributes other than
    metadata are lost.
    """
    @wraps(fn)
    def wrapper(sequence):
        return fn(sequence[:-1]) + sequence[-1]

    if isinstance(fn, ContextfulGenerator):

        @wraps(fn.learn_context)
        def learn_context(sequences: Iterable[SeqLike]):
            sequences = (sequence[:-1] for sequence in sequences)
            return fn.learn_context(sequences)
        wrapper.learn_context = learn_context  # type: ignore

    return wrapper


def keepsterm[T: SeqLike](fn: DecoyGenerator[T]) -> DecoyGenerator[T]:
    """Decorator that transform a DecoyGenerator into a new DecoyGenerator that
    doesn't alter the terminal aas.

    The `sequence` value is passed directly without the aminoacids that should
    be preserved. They are reinserted after the decorated function returns.

    :py:class:`strategies.ContextfulGenerator` objects preserve the
    :py:meth:`strategies.ContextfulGenerator.learn_context` method, but will
    receive each sequence without its terminal aas.

    Parameters
    ----------
    fn
        A :py:type:`strategies.DecoyGenerator`.

    Returns
    -------
    A version of `fn` that doesn't alter the terminal aminoacids of the target
    protein. If `fn` is a :py:class:`strategies.ContextfulGenerator`, the
    returned function will also be.

    Notes
    -----
    This function returns a new closure, meaning attributes other than
    metadata are lost.
    """
    def wrapper(sequence):
        return sequence[0] + fn(sequence[1:-1]) + sequence[-1]

    if isinstance(fn, ContextfulGenerator):

        @wraps(fn.learn_context)
        def learn_context(sequences: Iterable[SeqLike]):
            sequences = (sequence[1:-1] for sequence in sequences)
            return fn.learn_context(sequences)
        wrapper.learn_context = learn_context  # type: ignore

    return wrapper

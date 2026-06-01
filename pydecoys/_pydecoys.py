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

"""Internal file to specify the logic of the module's main API."""

from __future__ import annotations

from contextlib import contextmanager, nullcontext
import io
import os
import typing as t

try:
    _HAS_BIO = True
    import Bio  # noqa: F401
    from Bio.SeqRecord import SeqRecord
except ImportError as e:
    _HAS_BIO = False
    from warnings import warn
    warn(f"Module 'Biopython' not found: {str(e)}")

from pydecoys import strategies
from pydecoys.strategies import SeqLike
from pydecoys._builtins import decoy_strategy


type _PathOrIO = str | os.PathLike[str] | t.TextIO
type _Strategy = str | strategies.DecoyGenerator


def from_fasta(
    input: _PathOrIO,
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Iterable[SeqRecord]:
    """Apply a decoy generator to entries from a fasta file. Requires
    `Biopython`.

    Parameters
    ----------
    input
        Path or handle to a fasta file.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.

    Yields
    ------
    A decoy version of the next SeqRecord in the file.

    Notes
    -----
    This function calls :func:`from_SeqRecords` internally, so
    ``pydecoys.from_fasta(input, strategy)`` is equivalent to:

    >>> from Bio import SeqIO
    >>> targets = SeqIO.parse(input, format='fasta')  # doctest: +SKIP
    >>> decoys = from_SeqRecords(targets, strategy)   # doctest: +SKIP
    """

    if not _HAS_BIO:
        raise ImportError(
            "Module 'Biopython' necessary for IO operations. "
            "Install it with: 'pip install biopython'."
        )

    from Bio import SeqIO
    targets = SeqIO.parse(input, format='fasta')
    yield from from_SeqRecords(targets, strategy, decoy_tag, prefix)


def to_fasta(
    input: _PathOrIO | t.Iterable[SeqRecord] | SeqRecord,
    output: _PathOrIO,
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
    concat: bool = False
) -> int:
    """Apply a decoy generator to a set of sequences and write them to a file.
    Requires `Biopython`.

    Parameters
    ----------
    input
        A list (or iterator) of `SeqRecord` objects, a single `SeqRecord`, or
        the path or handle to a fasta file.
    output
        Path or handle to a fasta file.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.
    concat
        If `True`, output fasta will have both target and decoy entries. All
        targets will be first, then all decoys. Defaults to `False`.'

    Returns
    -------
    The number of entries written (as an integer).

    Notes
    -----
    While `_PathOrIO` is defined as
    ``str | os.PathLike[str] | typing.TextIO``, the runtime ``isinstance``
    check is done against ``io.TextIOBase``.

    This function calls :func:`from_SeqRecords` internally, so
    ``count = pydecoys.to_fasta(targets, output, strategy)`` where `targets`
    is an ``Iterable[SeqRecord]`` is equivalent to:

    >>> decoys = from_SeqRecords(targets, strategy)          # doctest: +SKIP
    >>> count = SeqIO.write(decoys, output, format='fasta')  # doctest: +SKIP

    If `True`, the `concat` flag will cause all sequences to be loaded
    into a list prior to decoy generation to avoid parsing sequences twice.
    """

    if not _HAS_BIO:
        raise ImportError(
            "Module 'Biopython' necessary for IO operations. "
            "Install it with: 'pip install biopython'."
        )

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    if isinstance(input, str | os.PathLike | io.TextIOBase):
        sequences = SeqIO.parse(input, 'fasta')
    elif isinstance(input, SeqRecord):
        sequences = [input]
    else:
        sequences = input

    if concat:
        from itertools import chain
        sequences = list(sequences)
        decoys = from_SeqRecords(sequences, strategy, decoy_tag, prefix)  # type: ignore
        write = chain(sequences, decoys)
    else:
        write = from_SeqRecords(sequences, strategy, decoy_tag, prefix)  # type: ignore

    return SeqIO.write(write, output, format='fasta')  # type: ignore


def from_SeqRecords(
    sequences: t.Iterable[SeqRecord] | SeqRecord,
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[SeqRecord, None, None]:
    """Apply a decoy generation strategy to a set of sequences. Requires
    `Biopython`.

    Parameters
    ----------
    sequences
        A list (or iterator) of `SeqRecord` objects, or a single `SeqRecord`.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.

    Yields
    ------
    A decoy version of the next SeqRecord in `sequences`.

    Examples
    --------
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> seqs = [
    ...     SeqRecord(Seq('DNIDYKAVYR'), 'seq1'),
    ...     SeqRecord(Seq('QSYMCTVTHP'), 'seq2'),
    ...     SeqRecord(Seq('CQWSLTEELL'), 'seq3'),
    ... ]
    >>> for decoy in from_SeqRecords(seqs, 'reverse'):
    ...     print(f'{decoy.id}: {decoy.seq}')
    decoy_seq1: RYVAKYDIND
    decoy_seq2: PHTVTCMYSQ
    decoy_seq3: LLEETLSWQC

    The `sequence` can be a single value:

    >>> from types import GeneratorType
    >>> seq = SeqRecord(Seq('DNIDYKAVYR'), 'seq1')
    >>> decoy = from_SeqRecords(seq, 'reverse')
    >>> isinstance(decoy, GeneratorType)  # Still returns a Generator
    True

    Notes
    -----
    This function will automatically pass the sequences as context to an unset
    :class:`strategies.ContextfulGenerator` strategy and reset it afterwards.
    """

    if not _HAS_BIO:
        raise ImportError(
            "Module 'Biopython' necessary for SeqRecord operations. "
            "Install it with: 'pip install biopython'."
        )

    from pydecoys import _bio
    sequences = _bio.iter_SeqRecord(sequences)
    tuples = (_bio.SeqRecord_to_tuple(record) for record in sequences)
    decoys = from_tuples(tuples, strategy, decoy_tag, prefix)
    records = (_bio.tuple_to_SeqRecord(decoy) for decoy in decoys)
    yield from records


def from_tuples[T: SeqLike](
    sequences: t.Iterable[tuple[str, T]],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[tuple[str, T], None, None]:
    """Apply a decoy generation strategy to a set of tuples.

    Differently from other functions in this module, `from_tuples` cannot
    accept a single tuple. Be sure to always pass an iterable of tuples.
    If you'd like, you can use :func:`tuple_as_decoy` instead.

    Parameters
    ----------
    sequences
        A list (or iterator) of `tuple` objects. The first item should be the
        seqid, and the second item should be the sequence.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.

    Yields
    ------
    A decoy version of the next tuple in `sequences`.

    Examples
    --------
    >>> seqs = [
    ...     ('seq1', 'DNIDYKAVYR'),
    ...     ('seq2', 'QSYMCTVTHP'),
    ...     ('seq3', 'CQWSLTEELL'),
    ... ]
    >>> for decoy in from_tuples(seqs, 'reverse'):
    ...     print(f'{decoy[0]}: {decoy[1]}')
    decoy_seq1: RYVAKYDIND
    decoy_seq2: PHTVTCMYSQ
    decoy_seq3: LLEETLSWQC

    Type is preserved:

    >>> from Bio.Seq import Seq
    >>> seqs = [
    ...     ('seq1', Seq('DNIDYKAVYR')),
    ...     ('seq2', Seq('QSYMCTVTHP')),
    ...     ('seq3', Seq('CQWSLTEELL')),
    ... ]
    >>> for decoy in from_tuples(seqs, 'reverse'):
    ...     print(isinstance(decoy[1], Seq))
    True
    True
    True

    Notes
    -----
    This function will automatically pass the sequences as context to an unset
    :class:`strategies.ContextfulGenerator` strategy and reset it afterwards.
    """
    decoy_generator = _validate_strategy(strategy)

    _validate_tag(decoy_tag)

    if isinstance(decoy_generator, strategies.ContextfulGenerator):
        sequences = list(sequences)

    with _get_contextualized(sequences, lambda x: x[1], decoy_generator) as gen:
        for seq in sequences:
            if not seq[1]:
                raise ValueError(f"Seq not present for sequence '{seq[0]}'")
            yield (_build_id(seq[0], decoy_tag, prefix), gen(seq[1]))


def from_seqs[T: SeqLike](
    sequences: t.Iterable[T] | T,
    strategy: _Strategy,
) -> t.Generator[T, None, None]:
    """Apply a decoy generation strategy to a set of sequences.

    Parameters
    ----------
    sequences
        A list (or iterator) of :type:`SeqLike` objects, or a single
        :type:`SeqLike`.
    strategy:
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.

    Yields
    ------
    A decoy version of the next sequence in `sequences`.

    Examples
    --------
    >>> seqs = [
    ...     'DNIDYKAVYR',
    ...     'QSYMCTVTHP',
    ...     'CQWSLTEELL',
    ... ]
    >>> for decoy in from_seqs(seqs, 'reverse'):
    ...     print(decoy)
    RYVAKYDIND
    PHTVTCMYSQ
    LLEETLSWQC

    Type is preserved:

    >>> from Bio.Seq import Seq
    >>> seqs = [
    ...     Seq('DNIDYKAVYR'),
    ...     Seq('QSYMCTVTHP'),
    ...     Seq('CQWSLTEELL'),
    ... ]
    >>> for decoy in from_seqs(seqs, 'reverse'):
    ...     print(isinstance(decoy, Seq))
    True
    True
    True

    The `sequence` can be a single value:

    >>> from types import GeneratorType
    >>> decoy = from_seqs('DNIDYKAVYR', 'reverse')
    >>> isinstance(decoy, GeneratorType)  # Still returns a Generator
    True

    Notes
    -----
    This function will automatically pass the sequences as context to an unset
    :class:`strategies.ContextfulGenerator` strategy and reset it afterwards.
    """

    decoy_generator = _validate_strategy(strategy)

    if _HAS_BIO:
        from Bio.Seq import Seq, MutableSeq
        if isinstance(sequences, str | Seq | MutableSeq):
            sequences = [sequences]  # type: ignore
    else:
        if isinstance(sequences, str):
            sequences = [sequences]  # type: ignore

    if isinstance(decoy_generator, strategies.ContextfulGenerator):
        sequences = list(sequences)  # type: ignore

    with _get_contextualized(sequences, lambda x: x, decoy_generator) as gen:
        for i, seq in enumerate(sequences):
            if not seq:
                raise ValueError(f'Seq not present for sequence {i}')
            yield gen(seq)  # type: ignore


def SeqRecord_as_decoy(
    sequence: SeqRecord,
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> SeqRecord:
    """Get a decoy from a given `SeqRecord`.

    Parameters
    ----------
    sequence
        A single `SeqRecord`.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.

    Returns
    -------
    A decoy version of `sequence`.

    Examples
    --------
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> seq = SeqRecord(Seq('DNIDYKAVYR'), 'seq1')
    >>> decoy = SeqRecord_as_decoy(seq, 'reverse')
    >>> print(f'{decoy.id}: {decoy.seq}')
    decoy_seq1: RYVAKYDIND

    Notes
    -----
    This function won't give context to :class:`strategies.ContextfulGerenator`
    objects. You can get an already set :class:`strategies.ContextfulGerenator`
    with :func:`get_contextualized_strategy`.

    """
    if not _HAS_BIO:
        raise ImportError(
            "Module 'Biopython' necessary for SeqRecord operations. "
            "Install it with: 'pip install biopython'."
        )

    from pydecoys import _bio
    seq_tuple = _bio.SeqRecord_to_tuple(sequence)
    decoy = tuple_as_decoy(seq_tuple, strategy, decoy_tag, prefix)
    return _bio.tuple_to_SeqRecord(decoy)


def tuple_as_decoy[T: SeqLike](
    sequence: tuple[str, T],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> tuple[str, T]:
    """Get a decoy from a given `tuple`.

    Parameters
    ----------
    sequence
        A single `tuple`. The first item should be the seqid, and the second
        item should be the sequence.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.

    Returns
    -------
    A decoy version of `sequence`.

    Examples
    --------
    >>> seq = ('seq1', 'DNIDYKAVYR')
    >>> decoy = tuple_as_decoy(seq, 'reverse')
    >>> print(f'{decoy[0]}: {decoy[1]}')
    decoy_seq1: RYVAKYDIND

    The type is preserved:

    >>> from Bio.Seq import Seq
    >>> seq = ('seq1', Seq('DNIDYKAVYR'))
    >>> decoy = tuple_as_decoy(seq, 'reverse')
    >>> isinstance(decoy[1], Seq)
    True

    Notes
    -----
    This function won't give context to :class:`strategies.ContextfulGerenator`
    objects. You can get an already set :class:`strategies.ContextfulGerenator`
    with :func:`get_contextualized_strategy`.
    """

    decoy_generator = _validate_strategy(strategy)

    _validate_tag(decoy_tag)

    if not sequence[1]:
        raise ValueError(f"Seq not present for sequence '{sequence[0]}'")

    if (
        isinstance(decoy_generator, strategies.ContextfulGenerator)
        and not decoy_generator.is_set
    ):
        raise ValueError(
            f"Strategy '{strategy}' requires context. "
            f"Use: `strategy = get_contextualized_strategy(sequences, {strategy})`"
        )

    if not sequence[1]:
        raise ValueError('Seq not present (cannot be empty str)')

    return (_build_id(sequence[0], decoy_tag, prefix), decoy_generator(sequence[1]))


def seq_as_decoy[T: SeqLike](
    sequence: T,
    strategy: _Strategy,
) -> T:
    """Get a decoy from a given :type:`SeqLike`.

    Parameters
    ----------
    sequence
        A single :type:`SeqLike`.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
         :type:`strategies.DecoyGenerator` function.

    Returns
    -------
    A decoy version of `sequence`.

    Examples
    --------
    >>> seq_as_decoy('DNIDYKAVYR', 'reverse')
    'RYVAKYDIND'

    The type is preserved:

    >>> from Bio.Seq import Seq
    >>> seq_as_decoy(Seq('DNIDYKAVYR'), 'reverse')
    Seq('RYVAKYDIND')

    Notes
    -----
    This function won't give context to :class:`strategies.ContextfulGerenator`
    objects. You can get an already set :class:`strategies.ContextfulGerenator`
    with :func:`get_contextualized_strategy`.
    """

    decoy_generator = _validate_strategy(strategy)

    if (
        isinstance(decoy_generator, strategies.ContextfulGenerator)
        and not decoy_generator.is_set
    ):
        raise ValueError(
            f"Strategy '{strategy}' requires context. "
            f"Use: `strategy = get_contextualized_strategy(sequences, {strategy})`"
        )

    if not sequence:
        raise ValueError('Seq not present (cannot be an empty str)')

    return decoy_generator(sequence)


def register(
    strategy_key: str,
    strategy_fn: strategies.DecoyGenerator
) -> None:
    """Register a new decoy strategy.

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
    >>> register('randomseq', random_seq)
    >>> seq = 'DNIDYKAVYR'
    >>> seq_as_decoy(seq, 'randomseq')  # doctest: +SKIP
    'LLEETLSWQC'

    The strategy key must be a lowercase string:

    >>> register(5, random_seq)
    Traceback (most recent call last):
        ...
    TypeError: Need a string for the decoy strategy (lower case)
    >>> register('RANDOMSEQ', random_seq)
    Traceback (most recent call last):
        ...
    ValueError: Strategy key 'RANDOMSEQ' should be lower case

    It must not be already defined:

    >>> register('randomseq', random_seq)
    Traceback (most recent call last):
        ...
    ValueError: Strategy key 'randomseq' already defined
    """

    if not isinstance(strategy_key, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy_key:
        raise ValueError("Strategy required (lower case string)")
    if not strategy_key.islower():
        raise ValueError(f"Strategy key '{strategy_key}' should be lower case")

    if strategy_key in decoy_strategy:
        raise ValueError(f"Strategy key '{strategy_key}' already defined")

    if not callable(strategy_fn):
        raise TypeError("Strategy function must be a callable")

    decoy_strategy[strategy_key] = strategy_fn


def get_contextualized_strategy(
    sequences: t.Iterable[SeqLike],
    strategy_key: str
) -> strategies.ContextfulGenerator:
    """Return a bare :type:`strategies.ContextfulGeneatir` from a
    context-based strategy key, with added context.

    Parameters
    ----------
    strategy_key
        Lower case string identifying the decoy strategy.

    Returns
    -------
    The correspondent decoy strategy.
    """

    import copy

    if not isinstance(strategy_key, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")

    strategy = _validate_strategy(strategy_key)

    if not isinstance(strategy, strategies.ContextfulGenerator):
        raise ValueError(f"Strategy '{strategy_key}' is not contextful.")

    strategy = copy.deepcopy(strategy)
    strategy.learn_context(sequences)

    return strategy


def _validate_strategy(strategy: _Strategy) -> strategies.DecoyGenerator:
    if callable(strategy):
        return strategy

    if not isinstance(strategy, str):
        raise TypeError("Need a string or callable for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string or callable)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    decoy_generator = decoy_strategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    return decoy_generator


def _validate_tag(decoy_tag: t.Any):
    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")


def _get_contextualized[T: SeqLike, U](
    sequences: t.Iterable[U],
    extract: t.Callable[[U], T],
    decoy_generator: strategies.DecoyGenerator[T],
):
    """Utilitary function.

    Normalizes `decoy_generator` as a ContextManager to ensure proper setup
    and cleanup for :clas:`strategies.ContextfulGenerator` instances that
    aren't set-up while avoiding code repetition for other cases.

    `extract` is just a callable to convert whatever type `U` is into a
    SeqLike.
    """
    if (
        isinstance(decoy_generator, strategies.ContextfulGenerator)
        and not decoy_generator.is_set
    ):
        extracted = (extract(seq) for seq in sequences)
        return _contextualized(decoy_generator, extracted)
    return nullcontext(decoy_generator)


def _build_id(seq_id: str, decoy_tag: str, prefix: bool) -> str:
    return decoy_tag + seq_id if prefix else seq_id + decoy_tag


@contextmanager
def _contextualized(
    strategy: strategies.ContextfulGenerator,
    sequences: t.Iterable[SeqLike]
):
    strategy.learn_context(sequences)
    try:
        yield strategy
    finally:
        strategy.reset()

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

import io
import os
import typing as t

if t.TYPE_CHECKING:
    from Bio.Seq import Seq, MutableSeq
    from Bio.SeqRecord import SeqRecord

from . import strategies
from .strategies import SeqLike


type _PathOrIO = str | os.PathLike[str] | t.TextIO
type _Strategy = str | strategies.DecoyGenerator


_decoy_strategy: dict[str, strategies.DecoyGenerator] = {
    "reverse": strategies.reverse,
    "reverse-keepn": strategies.reverse_keep_n,
    "reverse-keepc": strategies.reverse_keep_c,
    "reverse-keepterm": strategies.reverse_keep_term,
    "shuffle": strategies.shuffle,
    "shuffle-keepn": strategies.shuffle_keep_n,
    "shuffle-keepc": strategies.shuffle_keep_c,
    "shuffle-keepterm": strategies.shuffle_keep_term,
    "reversepep-trypsin": strategies.reversepep_trypsin,
    "reversepep-stricttrypsin": strategies.reversepep_stricttrypsin,
    "reversepep-argc": strategies.reversepep_argc,
    "reversepep-aspn": strategies.reversepep_aspn,
    "reversepep-chymo": strategies.reversepep_chymo,
    "reversepep-gluc": strategies.reversepep_gluc,
    "reversepep-lysc": strategies.reversepep_lysc,
    "reversepep-lysn": strategies.reversepep_lysn,
    "reversepep-stricttrypsin-keepn": strategies.reversepep_stricttrypsin_keepn,  # noqa: E501
    "shufflepep-trypsin": strategies.shufflepep_trypsin,
    "shufflepep-stricttrypsin": strategies.shufflepep_stricttrypsin,
    "shufflepep-argc": strategies.shufflepep_argc,
    "shufflepep-aspn": strategies.shufflepep_aspn,
    "shufflepep-chymo": strategies.shufflepep_chymo,
    "shufflepep-gluc": strategies.shufflepep_gluc,
    "shufflepep-lysc": strategies.shufflepep_lysc,
    "shufflepep-lysn": strategies.shufflepep_lysn,
    "shufflepep-stricttrypsin-keepn": strategies.shufflepep_stricttrypsin_keepn,  # noqa: E501
}


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
        :class:`strategies.DecoyGenerator`.
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
    >>> targets = SeqIO.parse(input, format='fasta')
    >>> decoys = pydecoys.from_SeqRecords(targets, strategy)
    """

    from Bio import SeqIO
    targets = SeqIO.parse(input, format='fasta')
    yield from from_SeqRecords(targets, strategy, decoy_tag, prefix)


def to_fasta(
    input: _PathOrIO | t.Iterable[SeqRecord] | SeqRecord,
    output: _PathOrIO,
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True
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
        :class:`strategies.DecoyGenerator`.
    decoy_tag
        An optional tag that is appended to each input's id. Defaults to
        `'decoy_'`.
    prefix
        If `True`, `decoy_tag` is prefixed, otherwise it's suffixed. Defaults
        to `True`.

    Returns
    -------
    The number of decoys written (as an integer).

    Notes
    -----
    While `_PathOrIO` is defined as
    ``str | os.PathLike[str] | typing.TextIO``, the runtime `isinstance` check
    is done against ``io.TextIOBase``.

    This function calls :func:`from_SeqRecords` internally, so
    ``count = pydecoys.to_fasta(targets, output, strategy)`` where `targets`
    is an ``Iterable[SeqRecord]`` is equivalent to:

    >>> from Bio import SeqIO
    >>> decoys = pydecoys.from_SeqRecords(targets, strategy)
    >>> count = SeqIO.write(decoys, output, format='fasta')
    """

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    if isinstance(input, str | os.PathLike | io.TextIOBase):
        sequences = from_fasta(input, 'fasta')
    elif isinstance(input, SeqRecord):
        sequences = [input]
    else:
        sequences = input

    decoys = from_SeqRecords(sequences, strategy, decoy_tag, prefix)  # type: ignore
    return SeqIO.write(decoys, output, format='fasta')


def from_SeqRecords(
    sequences: t.Iterable[SeqRecord] | SeqRecord,
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[SeqRecord, None, None]:
    """Apply a decoy generation strategy to a set of sequences.
    Requires `Biopython`.

    Parameters
    ----------
    sequences
        A list (or iterator) of `SeqRecord` objects, or a single `SeqRecord`.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
        :class:`strategies.DecoyGenerator`.
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
    >>> seqs = [
    ...     SeqRecord('DNIDYKAVYR', 'seq1'),
    ...     SeqRecord('QSYMCTVTHP', 'seq2'),
    ...     SeqRecord('CQWSLTEELL', 'seq3'),
    ... ]
    >>> for decoy in from_SeqRecords(seqs, 'reverse'):
    ...     print(f'{decoy.id}: {decoy.seq}')
    decoy_seq1: Seq('RYVAKYDIND')
    decoy_seq2: Seq('PHTVTCMYSQ')
    decoy_seq3: Seq('LLEETLSWQC')
    """

    from . import _bio
    sequences = _bio.iter_SeqRecord(sequences)
    tuples = (_bio.SeqRecord_to_tuple(record) for record in sequences)
    decoys = from_tuples(tuples, strategy, decoy_tag, prefix)
    records = (_bio.tuple_to_SeqRecord(decoy) for decoy in decoys)
    yield from records


@t.overload
def from_tuples(
    sequences: t.Iterable[tuple[str, str]],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[tuple[str, str], None, None]:
    ...


@t.overload
def from_tuples(
    sequences: t.Iterable[tuple[str, Seq]],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[tuple[str, Seq], None, None]:
    ...


@t.overload
def from_tuples(
    sequences: t.Iterable[tuple[str, MutableSeq]],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[tuple[str, MutableSeq], None, None]:
    ...


def from_tuples(
    sequences: t.Iterable[tuple[str, SeqLike]],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> t.Generator[tuple[str, SeqLike], None, None]:
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
        :class:`strategies.DecoyGenerator`.
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
    decoy_seq1: Seq('RYVAKYDIND')
    decoy_seq2: Seq('PHTVTCMYSQ')
    decoy_seq3: Seq('LLEETLSWQC')
    """

    decoy_generator = _validate_strategy(strategy)

    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")

    if isinstance(decoy_generator, strategies.ContextfulGenerator):
        sequences = list(sequences)
        # We extract the protein sequences itself
        seqs_only = [s[1] for s in sequences]
        decoy_generator.learn_context(seqs_only)

    for i, sequence in enumerate(sequences):
        if not sequence[1]:
            raise ValueError(f"Seq not present for tuple {i}: '{sequence[0]}'")

        id = decoy_tag + sequence[0] if prefix else sequence[0] + decoy_tag
        seq = decoy_generator(sequence[1])
        yield (id, seq)


@t.overload
def from_seqs(
    sequences: t.Iterable[str] | str,
    strategy: _Strategy,
) -> t.Generator[str, None, None]:
    ...


@t.overload
def from_seqs(
    sequences: t.Iterable[Seq] | Seq,
    strategy: _Strategy,
) -> t.Generator[Seq, None, None]:
    ...


@t.overload
def from_seqs(
    sequences: t.Iterable[MutableSeq] | MutableSeq,
    strategy: _Strategy,
) -> t.Generator[MutableSeq, None, None]:
    ...


def from_seqs(
    sequences: t.Iterable[SeqLike] | SeqLike,
    strategy: _Strategy,
) -> t.Generator[SeqLike, None, None]:
    """Apply a decoy generation strategy to a set of sequences.

    Parameters
    ----------
    sequences
        A list (or iterator) of :obj:`SeqLike` objects, or a single
        :obj:`SeqLike`.
    strategy:
        Lower case string specifying the decoy strategy to be used, or a
        :class:`strategies.DecoyGenerator`.

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
    'RYVAKYDIND'
    'PHTVTCMYSQ'
    'LLEETLSWQC'
    """

    decoy_generator = _validate_strategy(strategy)

    try:
        from Bio.Seq import Seq, MutableSeq
        if isinstance(sequences, str | Seq | MutableSeq):
            sequences = [sequences]

    except ModuleNotFoundError:
        if isinstance(sequences, str):
            sequences = [sequences]

    # if isinstance(decoy_generator, strategies.StatefulGenerator):
    #     sequences = decoy_generator.learn_state(sequences)

    if isinstance(decoy_generator, strategies.ContextfulGenerator):
        sequences = list(sequences)
        seqs_only = [s[1] for s in sequences]
        decoy_generator.learn_context(seqs_only)

    for i, sequence in enumerate(sequences):
        if not sequence:
            raise ValueError(f"No seq present for item {i}")

        yield decoy_generator(sequence)


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
        :class:`strategies.DecoyGenerator`.
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
    >>> seq = SeqRecord('DNIDYKAVYR', 'seq1')
    >>> decoy = SeqRecord_as_decoy(seq, 'reverse')
    >>> print(f'{decoy.id}: {decoy.seq}')
    decoy_seq1: Seq('RYVAKYDIND')
    """

    decoy_generator = _validate_strategy(strategy)

    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")

    if sequence.seq is None:
        raise ValueError(f"Seq not present for SeqRecord '{sequence.id}'")

    id = sequence.id if sequence.id else ""
    id = decoy_tag + id if prefix else id + decoy_tag
    seq = decoy_generator(sequence.seq)
    return SeqRecord(seq, id, description="")


@t.overload
def tuple_as_decoy(
    sequence: tuple[str, str],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> tuple[str, str]:
    ...


@t.overload
def tuple_as_decoy(
    sequence: tuple[str, Seq],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> tuple[str, Seq]:
    ...


@t.overload
def tuple_as_decoy(
    sequence: tuple[str, MutableSeq],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> tuple[str, MutableSeq]:
    ...


def tuple_as_decoy(
    sequence: tuple[str, SeqLike],
    strategy: _Strategy,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> tuple[str, SeqLike]:
    """Get a decoy from a given `tuple`.

    Parameters
    ----------
    sequence
        A single `tuple`. The first item should be the seqid, and the second
        item should be the sequence.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
        :class:`strategies.DecoyGenerator`.
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
    decoy_seq1: 'RYVAKYDIND'
    """

    decoy_generator = _validate_strategy(strategy)

    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")

    if not sequence[1]:
        raise ValueError(f"Seq not present for '{sequence[0]}'")

    id = decoy_tag + sequence[0] if prefix else sequence[0] + decoy_tag
    seq = decoy_generator(sequence[1])
    return (id, seq)


@t.overload
def seq_as_decoy(
    sequence: str,
    strategy: _Strategy,
) -> str:
    ...


@t.overload
def seq_as_decoy(
    sequence: Seq,
    strategy: _Strategy,
) -> Seq:
    ...


@t.overload
def seq_as_decoy(
    sequence: MutableSeq,
    strategy: _Strategy,
) -> MutableSeq:
    ...


def seq_as_decoy(
    sequence: SeqLike,
    strategy: _Strategy,
) -> SeqLike:
    """Get a decoy from a given :obj:`SeqLike`.

    Parameters
    ----------
    sequence
        A single :obj:`SeqLike`.
    strategy
        Lower case string specifying the decoy strategy to be used, or a
        :class:`strategies.DecoyGenerator`.

    Returns
    -------
    A decoy version of `sequence`.

    Examples
    --------
    >>> seq_as_decoy('DNIDYKAVYR', 'reverse')
    'RYVAKYDIND'
    """

    decoy_generator = _validate_strategy(strategy)
    return decoy_generator(sequence[1])


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
        A function following the :class:`strategies.DecoyGenerator` protocol.

    Examples
    --------
    Given a `random_seq` function that takes a sequence and returns a new,
    unrelated sequence of same size:

    >>> register('randomseq', random_seq)
    >>> seq = 'DNIDYKAVYR'
    >>> seq_as_decoy(seq, 'randomseq')
    'LLEETLSWQC'
    """

    if not isinstance(strategy_key, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy_key:
        raise ValueError("Strategy required (lower case string)")
    if not strategy_key.islower():
        raise ValueError(f"Strategy key '{strategy_key}' should be lower case")

    if strategy_key in _decoy_strategy:
        raise ValueError(f"Strategy '{strategy_key}' already exists")

    _decoy_strategy[strategy_key] = strategy_fn


def _validate_strategy(strategy: _Strategy) -> strategies.DecoyGenerator:
    if callable(strategy):
        return strategy

    if not isinstance(strategy, str):
        raise TypeError("Need a string or callable for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string or callable)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    decoy_generator = _decoy_strategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    return decoy_generator

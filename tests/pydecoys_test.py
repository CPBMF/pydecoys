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

# Most of this package's functionality is covered either by example docstrings
# (with doctest) or on an underlayer through `_builtins_test.py`

import importlib
import itertools

import pytest
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

from pydecoys._builtins import decoy_strategy
import pydecoys
from pydecoys._pydecoys import _validate_strategy, _validate_tag


KEYS = decoy_strategy.keys()

ID = 'test_sequence'
STR_SEQ = 'DNIDYKAVYR'
BIO_SEQ = Seq('DNIDYKAVYR')
MUT_SEQ = MutableSeq('DNIDYKAVYR')
STR_TUP = (ID, STR_SEQ)
BIO_TUP = (ID, BIO_SEQ)
MUT_TUP = (ID, MUT_SEQ)
SEQ_RECORD = SeqRecord(BIO_SEQ, id=ID)

DEC_ID = 'decoy_' + ID
DEC_STR = STR_SEQ[::-1]
DEC_BIO = Seq(STR_SEQ[::-1])
DEC_MUT = MutableSeq(STR_SEQ[::-1])
DEC_STR_TUP = (DEC_ID, DEC_STR)
DEC_BIO_TUP = (DEC_ID, DEC_BIO)
DEC_MUT_TUP = (DEC_ID, DEC_MUT)
DEC_RECORD = SeqRecord(DEC_BIO, id=DEC_ID)

PATH = 'tests/data/2026_01_ccp_crap.fasta'
CORRECT = 'tests/data/out/2026_01_ccp_crap_reverse.fasta'


@pytest.fixture
def fasta_path():
    return PATH


@pytest.fixture
def fasta_correct_path():
    return CORRECT


@pytest.fixture
def fasta_handle():
    return open(PATH, 'r')


@pytest.fixture
def fasta_correct_handle():
    return open(CORRECT, 'r')


@pytest.fixture
def tmp_file(tmp_path):
    return tmp_path / 'tmp_file.txt'


@pytest.fixture
def blank_handle(tmp_path):
    return open(tmp_path / 'blank.txt', 'w')


class DummyContextfulGenerator:
    def __init__(self):
        self.learned = False

    def learn_context(self, _):
        self.learned = True

    def __call__(self, _):
        pass


def dummy_strategy[T: pydecoys.strategies.SeqLike](_: T) -> T:
    raise NotImplementedError


def test_optional_Bio(missing_modules):
    with missing_modules('Bio'):
        importlib.reload(pydecoys)


@pytest.mark.parametrize('input', ['fasta_path', 'fasta_handle'])
def test_from_fasta(input, request):
    input = request.getfixturevalue(input)
    decoys = pydecoys.from_fasta(input, 'reverse')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


# Testing IO inputs
@pytest.mark.parametrize(
    ['input', 'concat'],
    itertools.product(
        ['fasta_path', 'fasta_handle'],
        [True, False]
    )
)
def test_to_fasta_ios(input, tmp_file, concat, request):
    input = request.getfixturevalue(input)
    pydecoys.to_fasta(PATH, tmp_file, 'reverse', concat=concat)
    decoys = SeqIO.parse(tmp_file, 'fasta')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    if concat:
        corrects = itertools.chain(SeqIO.parse(PATH, 'fasta'), corrects)
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


# No-IO inputs
@pytest.mark.parametrize(
    ['input', 'concat'],
    itertools.product(
        [SEQ_RECORD, [SEQ_RECORD]],
        [True, False]
    )
)
def test_to_fasta_not_IO(input, concat, tmp_file):
    pydecoys.to_fasta(input, tmp_file, 'reverse', concat=concat)
    decoys = SeqIO.parse(tmp_file, 'fasta')
    corrects = [DEC_RECORD]
    if concat:
        corrects = itertools.chain([SEQ_RECORD], corrects)
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


def test_from_SeqRecords():
    targets = SeqIO.parse(PATH, 'fasta')
    decoys = pydecoys.from_SeqRecords(targets, 'reverse')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


@pytest.mark.parametrize(
    ['input', 'corrects'],
    [
        ([STR_TUP], [DEC_STR_TUP]),
        ([BIO_TUP], [DEC_BIO_TUP]),
        ([MUT_TUP], [DEC_MUT_TUP])
    ]
)
def test_from_tuples(input, corrects):
    decoys = pydecoys.from_tuples(input, 'reverse')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert type(decoy[1]) is type(correct[1])
        assert decoy[0] == correct[0]
        assert decoy[1] == correct[1]


def test_from_tuples_contextful():
    contextful = DummyContextfulGenerator()
    sequences = [STR_TUP]
    list(pydecoys.from_tuples(sequences, contextful))
    assert contextful.learned


def test_from_tuples_empty():
    sequences = [STR_TUP, ('test_seq', '')]
    with pytest.raises(ValueError, match=r"Seq not present for sequence '.+'"):
        list(pydecoys.from_tuples(sequences, 'reverse'))


@pytest.mark.parametrize(
    ['input', 'corrects'],
    [
        (STR_SEQ,   [DEC_STR]),
        ([STR_SEQ], [DEC_STR]),
        (BIO_SEQ,   [DEC_BIO]),
        ([BIO_SEQ], [DEC_BIO]),
        (MUT_SEQ,   [DEC_MUT]),
        ([MUT_SEQ], [DEC_MUT])
    ]
)
def test_from_seqs(input, corrects):
    decoys = pydecoys.from_seqs(input, 'reverse')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert type(decoy) is type(correct)
        assert decoy == correct


@pytest.mark.parametrize(
    ['input', 'corrects'],
    [
        (STR_SEQ, [DEC_STR]),
        ([STR_SEQ], [DEC_STR])
    ]
)
def test_from_seqs_no_Bio(input, corrects, missing_modules):
    with missing_modules('Bio'):
        decoys = pydecoys.from_seqs(input, 'reverse')
        for decoy, correct in zip(decoys, corrects, strict=True):
            assert type(decoy) is type(correct)
            assert decoy == correct


def test_from_seqs_contextful():
    contextful = DummyContextfulGenerator()
    sequences = [STR_SEQ]
    list(pydecoys.from_seqs(sequences, contextful))
    assert contextful.learned


def test_from_seqs_empty():
    sequences = [STR_SEQ, '']
    with pytest.raises(ValueError, match=r'Seq not present for sequence \d+'):
        list(pydecoys.from_seqs(sequences, 'reverse'))


def test_SeqRecord_as_decoy():
    record = pydecoys.SeqRecord_as_decoy(SEQ_RECORD, 'reverse')
    assert record.id == DEC_RECORD.id
    assert record.seq == DEC_RECORD.seq


@pytest.mark.parametrize(
    ['input', 'correct'],
    [
        (STR_TUP, DEC_STR_TUP),
        (BIO_TUP, DEC_BIO_TUP),
        (MUT_TUP, DEC_MUT_TUP)
    ]
)
def test_tuple_as_decoy(input, correct):
    decoy = pydecoys.tuple_as_decoy(input, 'reverse')
    assert type(decoy) is type(correct)
    assert decoy[0] == correct[0]
    assert decoy[1] == correct[1]


def test_tuple_as_decoy_empty():
    with pytest.raises(ValueError, match=r"Seq not present for sequence '.+'"):
        pydecoys.tuple_as_decoy((ID, ''), 'reverse')


@pytest.mark.parametrize(
    ['input', 'correct'],
    [
        (STR_SEQ, DEC_STR),
        (BIO_SEQ, DEC_BIO),
        (MUT_SEQ, DEC_MUT)
    ]
)
def test_seq_as_decoy(input, correct):
    decoy = pydecoys.seq_as_decoy(input, 'reverse')
    assert type(decoy) is type(correct)
    assert decoy == correct


def test_seq_as_decoy_empty():
    with pytest.raises(ValueError, match=r"Seq not present \(cannot be an empty str\)"):
        pydecoys.seq_as_decoy('', 'reverse')


def test_validate_decoy_tag_bad_type():
    with pytest.raises(TypeError):
        _validate_tag(1)  # type: ignore


def test_register(monkeypatch):
    monkeypatch.setattr('pydecoys._pydecoys.decoy_strategy', {})
    pydecoys.register('dummy', dummy_strategy)
    strategy = pydecoys._pydecoys.decoy_strategy['dummy']  # type: ignore
    assert strategy is dummy_strategy


@pytest.mark.parametrize(
    ['key', 'fn'],
    [('dummy', 1), (1, dummy_strategy), (1, 1)]
)
def test_register_bad_type(key, fn):
    with pytest.raises(TypeError):
        pydecoys.register(key, fn)


def test_register_empty_str():
    matches = r"Strategy required \(lower case string\)"
    with pytest.raises(ValueError, match=matches):
        pydecoys.register('', dummy_strategy)


def test_register_str_not_lower():
    matches = r"Strategy key '.+' should be lower case"
    with pytest.raises(ValueError, match=matches):
        pydecoys.register('DUMMY', dummy_strategy)


def test_register_already_defined():
    matches = r"Strategy key '.+' already defined"
    with pytest.raises(ValueError, match=matches):
        pydecoys.register('reverse', dummy_strategy)


@pytest.mark.parametrize('strategy', ['reverse', dummy_strategy])
def test_validate_strategy(strategy):
    decoy_strategy = _validate_strategy(strategy)
    assert callable(decoy_strategy)


def test_validate_strategy_bad_type():
    with pytest.raises(TypeError):
        _validate_strategy(1)  # type: ignore


def test_validate_strategy_empty():
    matches = r"Strategy required \(lower case string or callable\)"
    with pytest.raises(ValueError, match=matches):
        _validate_strategy('')


def test_validate_strategy_upper():
    matches = r"Strategy string '.+' should be lower case"
    with pytest.raises(ValueError, match=matches):
        _validate_strategy('REVERSE')


def test_validate_strategy_unknown():
    matches = r"Unknown strategy: '.+'"
    with pytest.raises(ValueError, match=matches):
        _validate_strategy('unknown')

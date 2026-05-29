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
from typing import NoReturn

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

DEC_ID_SUFFIX = ID + '_decoy'
DEC_ID_CUSTOM = 'rev_' + ID
DEC_ID_CUSTOM_SUFFIX = ID + '_rev'

PATH = 'tests/data/2026_01_ccp_crap.fasta'
CORRECT = 'tests/data/out/2026_01_ccp_crap_reverse.fasta'

COUNT = len(list(SeqIO.parse(PATH, 'fasta')))


@pytest.fixture
def fasta_path():
    return PATH


@pytest.fixture
def fasta_handle():
    with open(PATH, 'r') as f:
        yield f


@pytest.fixture
def fasta_correct_handle():
    with open(CORRECT, 'r') as f:
        yield f


@pytest.fixture
def tmp_file(tmp_path):
    return tmp_path / 'tmp_file.txt'


class DummyContextfulGenerator:
    """Dummy class to emulate a ContextfulGenerator."""
    def __init__(self):
        self.learned = False

    def learn_context(self, _):
        self.learned = True

    def __call__(self, sequence):
        return sequence


def dummy_strategy(_) -> NoReturn:
    """Dummy function to emulate a DecoyGenerator."""
    raise NotImplementedError


def assert_decoy_equal(decoy, correct):
    """Helper function to assert that two seq representations are equal, no
    matter the actual type.
    """
    if isinstance(decoy, SeqRecord):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq
    elif isinstance(decoy, tuple):
        assert type(decoy[1]) is type(correct[1])
        assert decoy[0] == correct[0]
        assert decoy[1] == correct[1]
    else:
        assert type(decoy) is type(correct)
        assert decoy == correct


# Ensure code works without Biopython
# ===================================

def test_optional_Bio(missing_modules):
    with missing_modules('Bio'):
        importlib.reload(pydecoys)


# This is the only function that actually has a custom path without Biopython
# apart from raising an exception
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
            assert_decoy_equal(decoy, correct)


# Test decoy tagging
# ==================

@pytest.mark.parametrize(
    ['fn', 'input', 'correct_id'],
    [
        (pydecoys.SeqRecord_as_decoy, SEQ_RECORD, DEC_ID_SUFFIX),
        (pydecoys.tuple_as_decoy,     STR_TUP,    DEC_ID_SUFFIX),
    ]
)
def test_suffix_tag(fn, input, correct_id):
    result = fn(input, 'reverse', decoy_tag='_decoy', prefix=False)
    result_id = result.id if hasattr(result, 'id') else result[0]
    assert result_id == correct_id


@pytest.mark.parametrize(
    ['fn', 'input', 'correct_id'],
    [
        (pydecoys.SeqRecord_as_decoy, SEQ_RECORD, DEC_ID_CUSTOM),
        (pydecoys.tuple_as_decoy,     STR_TUP,    DEC_ID_CUSTOM),
    ]
)
def test_custom_tag_prefix(fn, input, correct_id):
    result = fn(input, 'reverse', decoy_tag='rev_')
    result_id = result.id if hasattr(result, 'id') else result[0]
    assert result_id == correct_id


@pytest.mark.parametrize(
    ['fn', 'input', 'correct_id'],
    [
        (pydecoys.SeqRecord_as_decoy, SEQ_RECORD, DEC_ID_CUSTOM_SUFFIX),
        (pydecoys.tuple_as_decoy,     STR_TUP,    DEC_ID_CUSTOM_SUFFIX),
    ]
)
def test_custom_tag_suffix(fn, input, correct_id):
    result = fn(input, 'reverse', decoy_tag='_rev', prefix=False)
    result_id = result.id if hasattr(result, 'id') else result[0]
    assert result_id == correct_id


# Testing IO functions
# ====================

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
def test_to_fasta_path(input, concat, tmp_file, request):
    input = request.getfixturevalue(input)

    orig = SeqIO.parse(PATH, 'fasta')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    corrects = _if_concat(orig, corrects, concat)

    count = pydecoys.to_fasta(input, tmp_file, 'reverse', concat=concat)
    decoys = SeqIO.parse(tmp_file, 'fasta')
    assert count == (COUNT*2 if concat else COUNT)
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert_decoy_equal(decoy, correct)


# Non-IO inputs
@pytest.mark.parametrize(
    ['input', 'concat'],
    itertools.product(
        [SEQ_RECORD, [SEQ_RECORD]],
        [True, False]
    )
)
def test_to_fasta_not_IO(input, concat, tmp_file):
    corrects = [DEC_RECORD]
    corrects = _if_concat([SEQ_RECORD], corrects, concat)

    count = pydecoys.to_fasta(input, tmp_file, 'reverse', concat=concat)

    assert count == (2 if concat else 1)
    decoys = SeqIO.parse(tmp_file, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert_decoy_equal(decoy, correct)


# Testing outputs
@pytest.mark.parametrize('concat', [True, False])
def test_to_fasta_output(fasta_path, tmp_file, concat):

    orig = SeqIO.parse(PATH, 'fasta')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    corrects = list(_if_concat(orig, corrects, concat))

    count = pydecoys.to_fasta(fasta_path, tmp_file, 'reverse', concat=concat)
    decoys = SeqIO.parse(tmp_file, 'fasta')
    assert count == (COUNT*2 if concat else COUNT)
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert_decoy_equal(decoy, correct)

    with open(tmp_file, 'w') as handle:
        count = pydecoys.to_fasta(fasta_path, handle, 'reverse', concat=concat)
    decoys = SeqIO.parse(tmp_file, 'fasta')
    assert count == (COUNT*2 if concat else COUNT)
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert_decoy_equal(decoy, correct)


def _if_concat(orig, corrects, concat):
    if concat:
        corrects = itertools.chain(orig, corrects)
    return corrects


# Testing the generator functions
# ===============================

@pytest.mark.parametrize(
    ['fn', 'inputs', 'corrects'],
    [
        (pydecoys.from_SeqRecords, [SEQ_RECORD],  [DEC_RECORD]),
        (pydecoys.from_SeqRecords, SEQ_RECORD,    [DEC_RECORD]),
        (pydecoys.from_tuples,     [STR_TUP],     [DEC_STR_TUP]),
        (pydecoys.from_tuples,     [BIO_TUP],     [DEC_BIO_TUP]),
        (pydecoys.from_tuples,     [MUT_TUP],     [DEC_MUT_TUP]),
        (pydecoys.from_seqs,       [STR_SEQ],     [DEC_STR]),
        (pydecoys.from_seqs,       [BIO_SEQ],     [DEC_BIO]),
        (pydecoys.from_seqs,       [MUT_SEQ],     [DEC_MUT]),
        (pydecoys.from_seqs,       STR_SEQ,       [DEC_STR]),
        (pydecoys.from_seqs,       BIO_SEQ,       [DEC_BIO]),
        (pydecoys.from_seqs,       MUT_SEQ,       [DEC_MUT]),
    ]
)
def test_from_batch(fn, inputs, corrects):
    decoys = fn(inputs, 'reverse')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert_decoy_equal(decoy, correct)


# Testing the single data functions
# =================================

@pytest.mark.parametrize(
    ['fn', 'input', 'correct'],
    [
        (pydecoys.SeqRecord_as_decoy, SEQ_RECORD, DEC_RECORD),
        (pydecoys.tuple_as_decoy,     STR_TUP,    DEC_STR_TUP),
        (pydecoys.tuple_as_decoy,     BIO_TUP,    DEC_BIO_TUP),
        (pydecoys.tuple_as_decoy,     MUT_TUP,    DEC_MUT_TUP),
        (pydecoys.seq_as_decoy,       STR_SEQ,    DEC_STR),
        (pydecoys.seq_as_decoy,       BIO_SEQ,    DEC_BIO),
        (pydecoys.seq_as_decoy,       MUT_SEQ,    DEC_MUT),
    ]
)
def test_as_decoy(fn, input, correct):
    decoy = fn(input, 'reverse')
    assert_decoy_equal(decoy, correct)


# Testing if ContextfulGenerators actually learn
# ==============================================

@pytest.mark.parametrize(
    ['fn', 'input'],
    [
        (pydecoys.from_SeqRecords, [SEQ_RECORD]),
        (pydecoys.from_tuples,     [STR_TUP]),
        (pydecoys.from_seqs,       [STR_SEQ]),
    ]
)
def test_contextful_learns(fn, input):
    contextful = DummyContextfulGenerator()
    list(fn(input, contextful))
    assert contextful.learned


# Empty seqs
# ==========

def test_from_tuples_empty():
    sequences = [STR_TUP, ('test_seq', '')]
    with pytest.raises(ValueError, match=r"Seq not present for sequence '.+'"):
        list(pydecoys.from_tuples(sequences, 'reverse'))


def test_from_seqs_empty():
    sequences = [STR_SEQ, '']
    with pytest.raises(ValueError, match=r'Seq not present for sequence [0-9]+'):
        list(pydecoys.from_seqs(sequences, 'reverse'))


def test_tuple_as_decoy_empty():
    with pytest.raises(ValueError, match=r"Seq not present for sequence '.+'"):
        pydecoys.tuple_as_decoy((ID, ''), 'reverse')


def test_seq_as_decoy_empty():
    with pytest.raises(ValueError, match=r"Seq not present \(cannot be an empty str\)"):
        pydecoys.seq_as_decoy('', 'reverse')


def test_validate_decoy_tag_bad_type():
    with pytest.raises(TypeError):
        _validate_tag(1)  # type: ignore


# Register function
# =================

def test_register(monkeypatch):
    monkeypatch.setattr('pydecoys._pydecoys.decoy_strategy', {})
    pydecoys.register('dummy', dummy_strategy)
    strategy = pydecoys._pydecoys.decoy_strategy['dummy']  # type: ignore
    assert strategy is dummy_strategy


@pytest.mark.parametrize(
    ['key', 'fn'],
    [
        ('dummy', 1),         # Bad fn
        (1, dummy_strategy),  # Bad key
        (1, 1),               # Both bad
    ],
    ids=['bad_fn', 'bad_key', 'both_bad']
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


# Ensuring strategy keys are correctly validated
# ==============================================

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

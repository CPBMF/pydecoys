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

import pytest
from Bio import SeqIO

from pydecoys._builtins import decoy_strategy
import pydecoys
from pydecoys._pydecoys import _validate_strategy


KEYS = decoy_strategy.keys()
PATH = 'tests/data/2026_01_ccp_crap.fasta'
CORRECT = 'tests/data/out/2026_01_ccp_crap_reverse.fasta'


@pytest.fixture
def fasta_handle():
    return open('tests/data/2026_01_ccp_crap.fasta', 'r')


@pytest.fixture
def tmp_file(tmp_path):
    return tmp_path / 'tmp_file.txt'


@pytest.fixture
def blank_handle(tmp_path):
    return open(tmp_path / 'blank.txt', 'w')


@pytest.fixture
def SeqRecord_iter():
    return SeqIO.parse(PATH, 'fasta')


def test_from_fasta_path():
    decoys = pydecoys.from_fasta(PATH, 'reverse')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


def test_from_fasta_handle(fasta_handle):
    decoys = pydecoys.from_fasta(fasta_handle, 'reverse')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


def test_to_fasta_path(tmp_file):
    pydecoys.to_fasta(PATH, tmp_file, 'reverse')
    decoys = SeqIO.parse(tmp_file, 'fasta')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


def test_to_fasta_handle(fasta_handle, tmp_file):
    pydecoys.to_fasta(fasta_handle, tmp_file, 'reverse')
    decoys = SeqIO.parse(tmp_file, 'fasta')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


def test_to_fasta_iter(SeqRecord_iter, tmp_file):
    pydecoys.to_fasta(SeqRecord_iter, tmp_file, 'reverse')
    decoys = SeqIO.parse(tmp_file, 'fasta')
    corrects = SeqIO.parse(CORRECT, 'fasta')
    for decoy, correct in zip(decoys, corrects, strict=True):
        assert decoy.id == correct.id
        assert decoy.seq == correct.seq


def test_validate_strategy_bad_type():
    with pytest.raises(TypeError):
        _validate_strategy(1)  # type: ignore


def test_validate_strategy_empty():
    with pytest.raises(ValueError):
        _validate_strategy('')


def test_validate_strategy_upper():
    with pytest.raises(ValueError, match=r"Strategy string '.+' should be lower case"):
        _validate_strategy('REVERSE')


def test_validate_strategy_unknown():
    with pytest.raises(ValueError, match=r"Unknown strategy: '.+'"):
        _validate_strategy('unknown')

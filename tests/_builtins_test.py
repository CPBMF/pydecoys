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

import itertools

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import pytest

from pydecoys import _builtins as _b
from pydecoys.strategies import DecoyGenerator, RAND, ContextfulGenerator, SeqLike


SEQ = 'DNIDYKAVYR'
TYPES: list[type] = [str, Seq, MutableSeq]
KEYS = _b.decoy_strategy.keys()
FUNCS = _b.decoy_strategy.values()


@pytest.fixture(autouse=True)
def rng():
    RAND.seed(10)


# We ensure each builtin returns the correct type
@pytest.mark.parametrize(['type', 'fn'], itertools.product(TYPES, FUNCS))
def test_bultin_types(type: type[SeqLike], fn: DecoyGenerator):
    if isinstance(fn, ContextfulGenerator):
        fn.learn_context([SEQ])
    seq = fn(type(SEQ))
    assert isinstance(seq, type)


@pytest.mark.parametrize('key', KEYS)
def test_builtins(key: str):
    filename = key.replace('-', '_')
    fn = _b.decoy_strategy[key]

    targets = SeqIO.parse('tests/data/2026_01_ccp_crap.fasta', 'fasta')

    if isinstance(fn, ContextfulGenerator):
        target_seqs = [target.seq for target in targets]
        fn.learn_context(target_seqs)

    decoys = (fn(record.seq) for record in targets)
    corrects = SeqIO.parse(f'tests/data/out/2026_01_ccp_crap_{filename}.fasta', 'fasta')

    for decoy, correct in zip(decoys, corrects):
        assert decoy == correct.seq

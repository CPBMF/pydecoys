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

from Bio import SeqIO
import pytest

from pydecoys.core import _get_contextualized
from pydecoys import _builtins as _b


KEYS = _b.decoy_strategy.keys()


@pytest.mark.parametrize('key', KEYS)
def test_builtins(key, root):
    filename = key.replace('-', '_')
    fn = _b.decoy_strategy[key]

    targets = list(SeqIO.parse(root / 'data/2026_01_ccp_crap.fasta', 'fasta'))
    output = root / f'data/out/2026_01_ccp_crap_{filename}.fasta'

    with _get_contextualized(targets, lambda x: x.seq, fn) as cxt:
        decoys = [cxt(record.seq) for record in targets]
        corrects = list(SeqIO.parse(output, 'fasta'))

        for decoy, correct in zip(decoys, corrects, strict=True):
            assert decoy == correct.seq

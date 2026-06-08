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

import filecmp

import pytest

from pydecoys.__main__ import main
from pydecoys._builtins import DECOY_STRATEGIES


KEYS = DECOY_STRATEGIES.keys()
PATH = 'tests/data/2026_01_ccp_crap.fasta'


@pytest.mark.parametrize('strategy', KEYS)
def test_main(strategy, tmp_path):
    tmp_file = tmp_path / 'decoys.fasta'
    correct = f'tests/data/out/2026_01_ccp_crap_{strategy.replace('-', '_')}.fasta'
    args = [PATH, '-o', str(tmp_file), '-s', strategy]
    main(args)
    assert filecmp.cmp(tmp_file, correct)

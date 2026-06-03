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
from typing import override

import pytest

from pydecoys import strategies as s


class DummyEnzymeSpecificGenerator(s.EnzymeSpecificGenerator):
    @override
    def __call__[T: s.SeqLike](self, sequence: T) -> T:
        raise NotImplementedError


class TestEnzymeSpecificGenerator:
    GOOD_CUT = ['K', 'R', 'F', 'KR', 'FY', 'FWL']
    GOOD_NOCUT = ['D', 'E', 'M', 'DE', 'MA', 'MGH', None]
    GOOD_SENSE = ['N', 'C']

    SHARED_CUT = ['K', 'KY', 'KRY']
    SHARED_NOCUT = ['K', 'KY', 'KRY']
    SHARED_NOCUT_N = ['K', 'KY', 'KRY']

    def test_not_str_cut(self):
        with pytest.raises(TypeError):
            DummyEnzymeSpecificGenerator(1)  # type: ignore

    def test_empty_cut(self):
        with pytest.raises(ValueError):
            DummyEnzymeSpecificGenerator("")

    def test_cut_not_aa(self):
        with pytest.raises(ValueError):
            DummyEnzymeSpecificGenerator("KRÇY")

    @pytest.mark.parametrize('param', ['nocut', 'nocut_n'])
    def test_not_optional_str_nocut(self, param):
        param_dict = {param: 1}
        with pytest.raises(TypeError):
            DummyEnzymeSpecificGenerator("KR", **param_dict)  # type: ignore

    @pytest.mark.parametrize('param', ['nocut', 'nocut_n'])
    def test_empty_nocut(self, param):
        param_dict = {param: ""}
        with pytest.raises(ValueError):
            DummyEnzymeSpecificGenerator("KR", **param_dict)  # type: ignore

    @pytest.mark.parametrize(
        ['cut', 'nocut', 'param'],
        itertools.product(SHARED_CUT, SHARED_NOCUT, ['nocut', 'nocut_n'])
    )
    def test_shared_cut_nocut(self, cut, nocut, param):
        param_dict = {param: nocut}
        with pytest.raises(ValueError):
            DummyEnzymeSpecificGenerator(cut, **param_dict)

    @pytest.mark.parametrize('param', ['nocut', 'nocut_n'])
    def test_nocut_not_aa(self, param):
        param_dict = {param: "KRÇY"}
        with pytest.raises(ValueError):
            DummyEnzymeSpecificGenerator("KR", **param_dict)  # type: ignore

    @pytest.mark.parametrize('sense', ['NC', 'B', '', 4])
    def test_bad_sense(self, sense):
        with pytest.raises(TypeError):
            DummyEnzymeSpecificGenerator("KR", sense=sense)

    @pytest.mark.parametrize(
        ['cut', 'nocut', 'sense'],
        itertools.product(GOOD_CUT, GOOD_NOCUT, GOOD_SENSE)
    )
    def test_good_init(self, cut, nocut, sense):
        DummyEnzymeSpecificGenerator(cut, nocut, sense)

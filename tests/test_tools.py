# This file is part of a program that make prediction of active
# protein phosphorylation sites using machine learning
# Copyright (C) 2018  Zoé Brunet
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import unittest
from utils.tools import find_pattern


class TestTools(unittest.TestCase):

    def test_find_pattern(self):
        pattern = r"coincoin"
        seq_with_pattern = "toto dit coincoin"
        seq_without_pattern = "toto dit meuh"
        result_with_pattern = find_pattern(pattern, seq_with_pattern)
        result_without_pattern = find_pattern(pattern, seq_without_pattern)
        assert result_with_pattern
        assert not result_without_pattern


if __name__ == '__main__':
    unittest.main()

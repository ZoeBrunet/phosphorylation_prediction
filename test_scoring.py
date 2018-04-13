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
from utils import scoring, parse_args


class TestUtils(unittest.TestCase):
    def test_scorint_zero_score(self):
        parser = parse_args(['X', 'example.fasta'])
        score = scoring(parser.pattern, parser.file)
        for s in score:
            self.assertEqual(s, 0)

    def test_scoring_normal_score(self):
        parser = parse_args(['T', 'example.fasta'])
        score = scoring(parser.pattern, parser.file)
        for s in score:
            assert(s >= 0)
            assert (s < 1)

    def test_scoring_max_score(self):
        parser = parse_args(['A', 'example.fasta'])
        score = scoring(parser.pattern, parser.file)
        for s in score:
            assert (s >= 0)
        self.assertEqual(score[1], 1)


if __name__ == '__main__':
    unittest.main()


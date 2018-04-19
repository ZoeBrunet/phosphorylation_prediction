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
from utils import parse_args, scoring, find_pattern


class TestUtils(unittest.TestCase):
    # The alignment of example.fasta gives :
    #
    # Example1: ----------AACCGTTCA
    # Example2: ----------AAAAAGAAA
    # Example3: ----------GACCG----
    # Example4: ----------AAAATGTTT
    # Example5: ----------CACCGATG-
    # Example6: AAAAAAAAAAAAAAAAAAA

    global example
    example = '~/dev/phosphorylation_prediction/data/example.fasta'

    def test_scoring_zero_score(self):
        parser = parse_args(['X', example])
        score = scoring(parser.pattern, parser.file, parser.max_window)
        for s in score:
            self.assertEqual(s, 0)

    def test_scoring_normal_score(self):
        parser = parse_args(['T', example])
        score = scoring(parser.pattern, parser.file, parser.max_window)
        for s in score:
            assert(s >= 0)
            assert (s < 1)
        self.assertAlmostEqual(score[16], 3/6)

    def test_scoring_max_score(self):
        parser = parse_args(['A', example])
        score = scoring(parser.pattern, parser.file, parser.max_window)
        for s in score:
            assert (s >= 0)
        self.assertAlmostEqual(score[11], 1)

    # According to the alignment score should be:
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(3 * 6) + 1/(3 * 6),
    # 2/(3 * 6) + 2/(3 * 6), 2/(3 * 6) + 2/(3 * 6) + 1/(3 * 6),
    # 1/(3 * 6) + 1/(3 * 6) + 2/(3 * 6), 1/(3 * 6), 0]
    def test_scoring_regexpr(self):
        parser = parse_args(['A.G', example])
        score = scoring(parser.pattern, parser.file, parser.max_window)
        for i in range(0, 13):
            self.assertEqual(score[i], 0)
        self.assertEqual(score[18], 0)
        self.assertAlmostEqual(score[13], 2 / 18)
        self.assertAlmostEqual(score[14], 4 / 18)
        self.assertAlmostEqual(score[15], 3 / 18)
        self.assertAlmostEqual(score[16], 2 / 18)
        self.assertAlmostEqual(score[17], 1 / 18)

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

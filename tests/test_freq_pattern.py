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
from utils.parser import common_parser
from utils.score import freq_of_pattern


class TestFreqPattern(unittest.TestCase):
    # The alignment of example.fasta gives :
    #
    # Example1: ----------AACCGTTCA
    # Example2: ----------AAAAAGAAA
    # Example3: ----------GACCG----
    # Example4: ----------AAAATGTTT
    # Example5: ----------CACCGATG-
    # Example6: AAAAAAAAAAAAAAAAAAA

    global example, window
    example = '/home/zoe/dev/phosphorylation_prediction/align/example_align.fasta'
    window =[0, 18]

    def test_scoring_zero_score(self):
        score = freq_of_pattern('X', window, example)
        for s in score:
            self.assertEqual(s, 0)

    def test_scoring_normal_score(self):
        score = freq_of_pattern('T', window, example)
        for s in score:
            assert(s >= 0)
            assert (s < 1)
        self.assertAlmostEqual(score[16], 3/6)

    def test_scoring_max_score(self):
        score = freq_of_pattern('A', window, example)
        for s in score:
            assert (s >= 0)
        self.assertAlmostEqual(score[11], 1)

    # According to the alignment score should be:
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(3 * 6) + 1/(3 * 6),
    # 2/(3 * 6) + 2/(3 * 6), 2/(3 * 6) + 2/(3 * 6) + 1/(3 * 6),
    # 1/(3 * 6) + 1/(3 * 6) + 2/(3 * 6), 1/(3 * 6), 0]
    def test_scoring_regexpr(self):
        score = freq_of_pattern('A.G', window, example)
        for i in range(0, 13):
            self.assertEqual(score[i], 0)
        self.assertEqual(score[18], 0)
        self.assertAlmostEqual(score[13], 2 / 18)
        self.assertAlmostEqual(score[14], 4 / 18)
        self.assertAlmostEqual(score[15], 3 / 18)
        self.assertAlmostEqual(score[16], 2 / 18)
        self.assertAlmostEqual(score[17], 1 / 18)


if __name__ == '__main__':
    unittest.main()

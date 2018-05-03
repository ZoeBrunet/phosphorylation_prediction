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
import os
from utils.score import *
import math

global example
my_path = os.path.abspath(os.path.dirname(__file__))
example = '%s/data/align/example_align.fasta' % my_path


class TestScore(unittest.TestCase):
    # The alignment of example.fasta gives :
    #
    # Example1: ----------AACCGTTCA
    # Example2: ----------AAAAAGAAA
    # Example3: ----------GACCG----
    # Example4: ----------AAAATGTTT
    # Example5: ----------CACCGATG-
    # Example6: AAAAAAAAAAAAAAAAAAA

    global example, window
    my_path = os.path.abspath(os.path.dirname(__file__))
    example = '%s/data/align/example_align.fasta' % my_path
    window = [0, 18]

    def test_scoring_zero_score(self):
        freq = get_freq_of_pattern('X', window, example)
        nb_align = freq["nb_align"]
        score = freq["score"]
        for s in score:
            self.assertEqual(s, 0)
        self.assertEqual(nb_align, 6)

    def test_scoring_normal_score(self):
        freq = get_freq_of_pattern('T', window, example)
        nb_align = freq["nb_align"]
        score = freq["score"]
        for s in score:
            assert (s >= 0)
            assert (s < 1)
        self.assertAlmostEqual(score[16], 3 / 6)
        self.assertEqual(nb_align, 6)

    def test_scoring_max_score(self):
        freq = get_freq_of_pattern('A', window, example)
        nb_align = freq["nb_align"]
        score = freq["score"]
        for s in score:
            assert (s >= 0)
        self.assertAlmostEqual(score[11], 1)
        self.assertEqual(nb_align, 6)

    # According to the alignment score should be:
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(3 * 6) + 1/(3 * 6),
    # 2/(3 * 6) + 2/(3 * 6), 2/(3 * 6) + 2/(3 * 6) + 1/(3 * 6),
    # 1/(3 * 6) + 1/(3 * 6) + 2/(3 * 6), 1/(3 * 6), 0]
    def test_scoring_regexpr(self):
        freq = get_freq_of_pattern('A.G', window, example)
        nb_align = freq["nb_align"]
        score = freq["score"]
        for i in range(0, 13):
            self.assertEqual(score[i], 0)
        self.assertEqual(score[18], 0)
        self.assertAlmostEqual(score[13], 2 / 18)
        self.assertAlmostEqual(score[14], 4 / 18)
        self.assertAlmostEqual(score[15], 3 / 18)
        self.assertAlmostEqual(score[16], 2 / 18)
        self.assertAlmostEqual(score[17], 1 / 18)
        self.assertEqual(nb_align, 6)

    def test_get_information_content(self):
        IC = get_information_content([10, 15], example)
        fa = [4, 6, 3, 3, 2]
        fc = [1, 0, 3, 3, 0]
        fg = [1, 0, 0, 0, 3]
        ft = [0, 0, 0, 0, 1]
        q = 0.05
        ICtheo = 0
        lamb = lambda n: n * math.log(n / q, 10) if n != 0 else 0
        for f in [fa, fc, fg, ft]:
            f[:] = [x / 6 for x in f]
        for na, nc, ng, nt in zip(fa, fc, fg, ft):
            for n in [na, nc, ng, nt]:
                ICtheo += lamb(n)
        self.assertAlmostEquals(ICtheo, IC)

    def test_get_shanon_entropy(self):
        pssm = get_pssm(summary_align=get_align_info(example))
        shanon_entropy = get_shanon_entropy([10, 15], pssm)
        fa = [4, 6, 3, 3, 2]
        fc = [1, 0, 3, 3, 0]
        fg = [1, 0, 0, 0, 3]
        ft = [0, 0, 0, 0, 1]
        lamb = lambda n: - n * math.log(n, 2) if n != 0 else 0
        for f in [fa, fc, fg, ft]:
            f[:] = [x / 6 for x in f]
        for na, nc, ng, nt, shanon in zip(fa, fc, fg, ft, shanon_entropy):
            shanon_expected = 0
            for n in [na, nc, ng, nt]:
                shanon_expected += lamb(n)
            self.assertAlmostEquals(shanon, shanon_expected)


if __name__ == '__main__':
    unittest.main()




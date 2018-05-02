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


class TestInformationContent(unittest.TestCase):

    def test_get_pssm(self):
        summary_align = get_align_info(example)
        pssm = get_pssm(summary_align)
        for row in pssm:
            self.assertTrue(row["A"] >= 1)
            self.assertTrue(sum(row.values()) > 0)
            self.assertTrue(sum(row.values()) < 7)

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


if __name__ == '__main__':
    unittest.main()




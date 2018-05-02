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
import pandas as pd
from utils.enrich_csv import enrich_csv


class TestEnrichCsv(unittest.TestCase):

    global my_path, csv1, csv2
    my_path = os.path.abspath(os.path.dirname(__file__))
    csv1 = '%s/data/csv/csv1.csv' % my_path
    csv2 = '%s/data/csv/csv2.csv' % my_path

    def test_enrich_csv(self):
        created_csv = enrich_csv(csv1, csv2)
        expected_csv = '%s/data/csv/expected_csv.csv' % my_path
        a = pd.read_csv(created_csv, sep=";")
        b = pd.read_csv(expected_csv, sep=";")
        self.assertTrue(a.equals(b))
        os.remove(created_csv)


if __name__ == '__main__':
    unittest.main()

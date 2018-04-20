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
from convert_id import *


class TestConvertID(unittest.TestCase):

    def test_convert_single_id(self):
        uniprotidlist = ["O08539"]
        data = listgeneID(uniprotidlist)["out"][0]
        self.assertEqual(data['_id'], '30948')

    def test_missmatch_id(self):
        uniprotidlist = ["t0t0", "tata_de_toto", "O08539"]
        data = listgeneID(uniprotidlist)["out"]
        self.assertFalse('_id' in data[0])
        self.assertTrue(data[1]['notfound'])
        self.assertEqual(data[2]['_id'], '30948')

    def test_missmatch_id(self):
        geneID = '30948'
        json = request_gene_id(geneID)
        self.assertTrue(json["data"][0] == "EOG090406M5")


if __name__ == '__main__':
    unittest.main()




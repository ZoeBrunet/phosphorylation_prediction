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
from utils.import_csv import *


class TestImportCSV(unittest.TestCase):

    def test_convert_single_id(self):
        mg = mygene.MyGeneInfo()
        uniprotID = "O08539"
        position = 304
        code = "S"
        sequence = "ceci est un test"
        gene = Gene(uniprotID, position, code, sequence)
        gene._set_geneID(mg, 1, 1)
        self.assertEqual(gene._get_geneID(), '30948')

    def test_missmatch_id(self):
        mg = mygene.MyGeneInfo()
        uniprotID = "tata_de_toto"
        position = 304
        code = "S"
        sequence = "ceci est un autre test"
        gene = Gene(uniprotID, position, code, sequence)
        gene._set_geneID(mg, 1, 1)
        self.assertTrue(gene._get_geneID() is None)


    def test_cluster_id(self):
        mg = mygene.MyGeneInfo()
        uniprotID = "O08539"
        position = 304
        code = "S"
        sequence = "ceci est un autre test"
        gene = Gene(uniprotID, position, code, sequence)
        gene._set_geneID(mg, 1, 1)
        gene._set_cluster()
        self.assertTrue(gene._get_cluster() == "EOG090406M5")


if __name__ == '__main__':
    unittest.main()




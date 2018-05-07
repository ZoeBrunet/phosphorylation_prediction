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
import pandas
from utils.import_csv import *


class TestImportCSV(unittest.TestCase):

    global example, window
    my_path = os.path.abspath(os.path.dirname(__file__))
    example = '%s/data/csv/sample.csv' % my_path

    def test_convert_single_id(self):
        df = pandas.DataFrame({'acc': ["O08539"]})
        mg = mygene.MyGeneInfo()
        index = create_index(df, mg)
        self.assertEqual(str(index.values[0][1]), '30948')

    def test_missmatch_id(self):
        df = pandas.DataFrame({'acc': ["tata_de_toto"]})
        mg = mygene.MyGeneInfo()
        index = create_index(df, mg)
        self.assertEqual(str(index.values[0][1]), "None")

    def test_cluster_id(self):
        df = pandas.DataFrame({'acc': ["O08539"]})
        mg = mygene.MyGeneInfo()
        index = create_index(df, mg)
        self.assertEqual(str(index.values[0][3]), "EOG090406M5")

    def test_neg_geneID_list(self):
        df = import_csv(example)
        mg = mygene.MyGeneInfo()
        gen_list = gen_uniprot_id_list_neg(df, mg)
        for gene in gen_list:
            self.assertTrue(gene._get_code() == "S")
            self.assertTrue(abs(gene._get_position() - 304) > 50)
            self.assertTrue(gene._get_uniprotID() == 'O08539')


if __name__ == '__main__':
    unittest.main()




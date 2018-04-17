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

import mygene
import requests
from print_info_phospho_elm import import_csv


def gen_uniprot_id_list(csv):
    df = import_csv(csv)
    return df['acc'].value_counts().keys().tolist()


def listgeneID(list):
    mg = mygene.MyGeneInfo()
    return mg.querymany(list,
                        scopes='symbol,accession',
                        fields='uniprot')


def uniprot2geneID(csv):
    uniprotidlist = gen_uniprot_id_list(csv)
    return listgeneID(uniprotidlist)

def call_odb_api(geneID):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1' % geneID
    response = requests.get(request)
    return response.json()

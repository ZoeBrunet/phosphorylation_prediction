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
import os
from print_info_phospho_elm import import_csv


class Gene:

    def __init__(self, uniprotID, geneID):
        self.uniprotID = uniprotID
        self.geneID = geneID
        self.cluster = []

    def _set_cluster(self, cluster):
        self.cluster = cluster

    def _get_geneID(self):
        return self.geneID


def gen_uniprot_id_list(csv):
    df = import_csv(csv)
    return df['acc'].value_counts().keys().tolist()


def listgeneID(list):
    mg = mygene.MyGeneInfo()
    return mg.querymany(list,
                        scopes='symbol,accession',
                        fields='uniprot',
                        returnall=True)


def uniprot2geneID(csv):
    uniprotidlist = gen_uniprot_id_list(csv)
    return listgeneID(uniprotidlist)


def request_cluster_id(clusterID, path):
    request = "'http://www.orthodb.org/fasta?id=%s'" % clusterID
    name = "%s.fasta" % clusterID
    path2fastas = "%s/fastas" % path
    create_fastas_folder = "touch %s" % path2fastas
    request_api = "curl %s -o %s/%s" % (request, path2fastas, name)
    final_request = "%s ; %s" % (create_fastas_folder, request_api)
    os.system(final_request)


def request_gene_id(geneID):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1&singlecopy=1' \
              % geneID
    response = requests.get(request)
    return response.json()


def print_trace(i, length, request):
    print("request the orthodb API %s %s/%s = %s"
          % (request, str(i + 1), str(length),
             str(round(((i + 1) / length) * 100, 2)) + "%"))


def import_ortholog(path, file_name):
    csv = "%s/%s" % (path, file_name)
    index = uniprot2geneID(csv)
    id_gene = []
    clusterlist = []
    for g in index["out"]:
        if '_id' in g:
            id_gene.append(Gene(g['query'], g['_id']))
        else:
            id_gene.append(Gene(g['query'], None))
    length = len(id_gene)
    for i, id in enumerate(id_gene):
        print_trace(i, length, "gene id")
        id_gene = id._get_geneID()
        if id_gene is not None:
            clusters = request_gene_id(id_gene)
            id._set_cluster(clusters)
            for cluster in clusters["data"]:
                if cluster not in clusterlist:
                    clusterlist.append(cluster)
    length_cluster = len(clusterlist)
    for k, clusterid in enumerate(clusterlist):
        print_trace(k, length_cluster, "cluster id")
        request_cluster_id(clusterid, path)
    return id_gene


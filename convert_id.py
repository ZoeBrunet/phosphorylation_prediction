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
import argparse
from utils import common_parse
from print_info_phospho_elm import import_csv


def parse_args_align(args):
    parser = argparse.ArgumentParser(description='Run align-file to create'
                                                 'dataset from phospho.ELM dump')
    parser.add_argument('path',
                        help='Where is your folder ?')
    common_parse(parser)
    return parser.parse_args(args)


class Gene:
    def __init__(self, uniprotID, position, code, sequence):
        self.uniprotID = uniprotID
        self.position = position
        self.code = code
        self.sequence = sequence
        self.geneID = None
        self.taxID = None
        self.cluster = None

    def _get_uniprotID(self):
        return self.uniprotID

    def _get_taxID(self):
        return self.taxID

    def _get_geneID(self):
        return self.geneID

    def _get_cluster(self):
        return self.cluster

    def _get_position(self):
        return self.position

    def _get_code(self):
        return self.code

    def _get_sequence(self):
        return self.sequence

    def _set_geneID(self, mg, i, length):
        self.geneID = None
        print_trace(i, length, "request the orthodb API for gene id")
        response = mg.query(self.uniprotID,
                            scope='symbol,accession',
                            fields='uniprot')["hits"]
        if len(response):
            self.geneID = response[0]["_id"]

    def _set_taxID(self, mg):
        if self.geneID is not None:
            self.taxID = mg.getgene(self.geneID)['taxid']

    def _set_cluster(self):
        request = request_gene_id(self.geneID)
        if len(request["data"]):
            self.cluster = request["data"][0]

    def set_info(self, mg, i, length_gene_list):
        self._set_geneID(mg, i, length_gene_list)
        self._set_taxID(mg)
        self._set_cluster()


def print_trace(i, length, request):
    print("%s %s/%s = %s"
          % (request, str(i + 1), str(length),
             str(round(((i + 1) / length) * 100, 2)) + "%"))


def gen_uniprot_id_list(csv, pattern):
    df = import_csv(csv)
    genelist = []
    for acc, position, code, sequence in zip(df["acc"],
                                             df["position"],
                                             df["code"],
                                             df["sequence"]):
        unique = True
        if str(code) not in str(pattern):
            unique = False
        if len(genelist):
            for gene in genelist:
                if (((gene._get_uniprotID() == acc
                        and gene._get_position() == position
                        and gene._get_sequence() == sequence))
                        or str(code) not in str(pattern)):
                    unique = False
                    break
        if unique:
            genelist.append(Gene(acc, position, code, sequence))
            print("Import %s from the csv file" % acc)
    return list(set(genelist))


def request_gene_id(geneID):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1' \
              '&singlecopy=1&limit=1' % geneID
    response = requests.get(request)
    return response.json()


def request_cluster_id(clusterID, path):
    name = "%s.fasta" % clusterID
    path2fastas = "%s/fastas" % path
    request_odb = "'http://www.orthodb.org/fasta?id=%s'" % clusterID
    request_api = "curl %s -o %s/%s" % (request_odb, path2fastas, name)
    create_fastas_folder = "if [ ! -d %s ] ; " \
                           "then mkdir %s; " \
                           "fi" % (path2fastas, path2fastas)
    final_request = "%s; " \
                    "cd %s; " \
                    "if [ ! -f %s ] ; " \
                    "then %s;" \
                    "fi" % (create_fastas_folder,
                            path2fastas,
                            name,
                            request_api)
    os.system(final_request)


def import_ortholog(path, file_name, pattern):
    mg = mygene.MyGeneInfo()
    csv = "%s/%s" % (path, file_name)
    gene_list = gen_uniprot_id_list(csv, pattern)
    length_gene_list = len(gene_list)
    for i, gene in enumerate(gene_list):
        gene.set_info(mg, i, length_gene_list)
        if gene._get_cluster() is not None:
            request_cluster_id(gene._get_cluster(), path)
    return gene_list

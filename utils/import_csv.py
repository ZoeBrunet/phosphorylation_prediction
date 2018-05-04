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
import pandas as pd
from utils.tools import print_trace, find_pattern


class Gene:
    def __init__(self, uniprotID, position, code, sequence, phosphorylation_site):
        self.uniprotID = uniprotID
        self.position = position
        self.code = code
        self.sequence = sequence
        self.geneID = None
        self.taxID = None
        self.cluster = None
        self.phosphorylation_site = phosphorylation_site

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

    def _get_phosphorylation_site(self):
        return self.phosphorylation_site

    def set_info(self, index):
        self.geneID = index[0][1]
        self.taxID = index[0][2]
        self.cluster = index[0][3]


def import_csv(csv):
    df = pd.read_csv(csv)
    # Convert data into category
    for cat in df.columns:
        if cat != "position":
            df[cat] = df[cat].astype('category')
    return df


def gen_uniprot_id_list_neg(df, pattern):
    genelist = []
    for acc, position, code, sequence in zip(df["acc"],
                                             df["position"],
                                             df["code"],
                                             df["sequence"]):
        unique = True
        if str(code) in str(pattern):
            for m in find_pattern(pattern, sequence):
                new_position = round((m.end() + m.start() - 1)/2 + 1)
                if abs(new_position - position) <= 50:
                    unique = False
                if len(genelist):
                    for gene in genelist:
                        if(gene._get_uniprotID() == acc
                                and abs(gene._get_position() - position) <= 50):
                            genelist.remove(gene)
                        if (((gene._get_uniprotID() == acc
                              and gene._get_position() == position
                              and gene._get_sequence() == sequence))):
                            unique = False
                        if (((gene._get_uniprotID() == acc
                              and gene._get_position() == new_position
                              and gene._get_sequence() == sequence))):
                            unique = False
                            break
                if unique:
                    genelist.append(Gene(acc, new_position, pattern, sequence, False))
                    print("Import %s sites from the csv file" % acc)
    return list(set(genelist))


def gen_uniprot_id_list(df, pattern):
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
            genelist.append(Gene(acc, position, code, sequence, True))
            print("Import %s sites from the csv file" % acc)
    return list(set(genelist))


def request_gene_id(geneID):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1' \
              '&singlecopy=1&limit=1' % geneID
    response = requests.get(request)
    return response.json()


def request_cluster_id(clusterID, path):
    name = "%s.fasta" % clusterID
    path2fastas = "%s/fastas" % path
    path2file = "%s/%s" % (path2fastas, name)
    if not os.path.exists(path2fastas):
        os.mkdir(path2fastas)
    if not os.path.exists(path2file):
        request_odb = "'http://www.orthodb.org/fasta?id=%s'" % clusterID
        request_api = "curl %s -o %s" % (request_odb, path2file)
        os.system(request_api)


def create_index(df, mg):
    uniprot_id = df['acc'].value_counts().keys().tolist()
    resp = mg.querymany(uniprot_id,scope='symbol,accession',
                            fields='uniprot')
    colonnes = ["uniprotID", "geneID", "taxID", "clusterID"]
    lignes = []
    length = len(resp)
    for i, r in enumerate(resp):
        geneID = None
        taxID = None
        clusterID = None
        uniprotID = r["query"]
        if "_id" in r:
            geneID = r["_id"]
            if 'taxid' in mg.getgene(geneID):
                taxID = mg.getgene(geneID)['taxid']
            request = request_gene_id(geneID)
            if len(request["data"]):
                clusterID = request["data"][0]
        ligne = (uniprotID, geneID, taxID, clusterID)
        lignes.append(ligne)
        print_trace(i, length, "create index for %s" % uniprotID)
    df = pd.DataFrame(data=lignes, columns=colonnes)
    return df


def import_ortholog(csv, pattern, phospho_sites):
    path = os.path.dirname(os.path.dirname(csv))
    mg = mygene.MyGeneInfo()
    df = import_csv(csv)
    index = create_index(df, mg)
    gene_list = gen_uniprot_id_list(df, pattern) if phospho_sites \
        else gen_uniprot_id_list_neg(df, pattern)
    length = len(gene_list)
    for i, gene in enumerate(gene_list):
        print_trace(i, length, "convert uniprotID into geneID")
        gene.set_info(index[index.uniprotID == gene._get_uniprotID()].values)
        if gene._get_cluster() is not None:
            request_cluster_id(gene._get_cluster(), path)
    return gene_list

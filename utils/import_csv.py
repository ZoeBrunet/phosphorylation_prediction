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
from biothings_client import get_client
from utils.tools import print_trace, find_pattern, is_metazoan


class Gene:
    def __init__(self, uniprotID, position, code, sequence, phosphorylation_site):
        self.uniprotID = uniprotID
        self.position = position
        self.code = code
        self.sequence = sequence
        self.geneID = None
        self.taxID = None
        self.cluster = None
        self.metazoan = None
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

    def _get_metazoan(self):
        return self.metazoan

    def _get_phosphorylation_site(self):
        return self.phosphorylation_site

    def _set_sequence(self, sequence):
        self.sequence = sequence

    def set_info(self, index):
        self.geneID = index[0][1]
        self.taxID = index[0][2]
        self.cluster = index[0][3]
        self.metazoan = index[0][4]


def import_csv(csv):
    df = pd.read_csv(csv)
    # Convert data into category
    for cat in df.columns:
        if cat != "position":
            df[cat] = df[cat].astype('category')
    return df


def gen_uniprot_id_list_neg(liste, pattern):
    genelist = []
    length = len(liste)
    for i, gene1 in enumerate(liste):
        sequence = gene1._get_sequence()
        position = gene1._get_position()
        acc1 = gene1._get_uniprotID()
        unique = True
        for m in find_pattern(pattern, sequence):
            new_position = round((m.end() + m.start() - 1)/2)
            if abs(new_position - position) <= 50:
                unique = False
            if len(genelist):
                for gene2 in genelist:
                    acc2 = gene2._get_uniprotID()
                    if(acc2 == acc1
                            and abs(gene2._get_position() - position) <= 50):
                        genelist.remove(gene2)
                    if (((acc2 == acc1
                          and gene2._get_position() == new_position
                          and gene2._get_sequence() == sequence))):
                        unique = False
                        break
            if unique:
                genelist.append(Gene(acc1, new_position, pattern, sequence, False))
                print_trace(i, length, "Import %s sites from the csv file for negative dataset" % acc1)
    return list(set(genelist))


def gen_uniprot_id_list(df, pattern):
    length = len(df)
    genelist = []
    for i, (acc, position, code, sequence) in enumerate(zip(df["acc"],
                                                            df["position"],
                                                            df["code"],
                                                            df["sequence"])):
        unique = True
        if str(code) not in str(pattern):
            unique = False
        if len(genelist):
            for gene in genelist:
                if (((gene._get_uniprotID() == acc
                      and gene._get_position() == position - 1
                      and gene._get_sequence() == sequence))
                        or str(code) not in str(pattern)):
                    unique = False
                    break
        if unique:
            genelist.append(Gene(acc, position - 1, code, sequence, True))
            print_trace(i, length, "Import %s sites from the csv file for positive dataset" % acc)
    return list(set(genelist))


def request_gene_id(geneID, s):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1' \
              '&singlecopy=1&limit=1' % geneID
    response = s.get(request)
    if response.status_code == 200:
        return response.json()
    else:
        print("status code for %s = %s" % (request, response.status_code))
        return None


def request_cluster_id(clusterID, path, s):
    name = "%s.fasta" % clusterID
    path2fastas = "%s/fastas" % path
    path2file = "%s/%s" % (path2fastas, name)
    if not os.path.exists(path2fastas):
        os.mkdir(path2fastas)
    if not os.path.exists(path2file):
        request_odb = 'http://www.orthodb.org/fasta?id=%s' % clusterID
        resp = s.get(request_odb)
        request_api = "curl %s -o %s" % (request_odb, path2file)
        os.system(request_api)


def create_index(list, mg):
    uniprot_id = []
    mt = get_client("taxon")
    for gene in list:
        acc = gene._get_uniprotID()
        if acc not in uniprot_id:
            uniprot_id.append(acc)
    resp = mg.querymany(uniprot_id, scope='symbol,accession',
                        fields='uniprot, taxid', species="all")
    colonnes = ["uniprotID", "geneID", "taxID", "clusterID", "metazoan"]
    lignes = []
    length = len(resp)
    with requests.Session() as s:
        for i, r in enumerate(resp):
            geneID = None
            taxID = None
            clusterID = None
            metazoan = None
            uniprotID = r["query"]
            if 'taxid' in r:
                taxID = r["taxid"]
                metazoan = is_metazoan(taxID, mt)
            if "_id" in r:
                geneID = r["_id"]
                request = request_gene_id(geneID, s)
                if request is not None:
                    if len(request["data"]):
                        clusterID = request["data"][0]
                ligne = (uniprotID, geneID, taxID, clusterID, metazoan)
                lignes.append(ligne)
            print_trace(i, length, "create index for %s" % uniprotID)
        df = pd.DataFrame(data=lignes, columns=colonnes)
    return df


def fill_gene(gene_list, index, path):
    length = len(gene_list)
    with requests.Session() as s:
        for i, gene in enumerate(gene_list):
            print_trace(i, length, "convert uniprotID into geneID")
            ind = index[index.uniprotID == gene._get_uniprotID()].values
            if len(ind):
                gene.set_info(ind)
            if gene._get_position():
                if gene._get_cluster() is not None:
                    request_cluster_id(gene._get_cluster(), path, s)


def import_ortholog(csv, pattern):
    path = os.path.dirname(os.path.dirname(csv))
    mg = mygene.MyGeneInfo()
    df = import_csv(csv)
    gene_list_pos = gen_uniprot_id_list(df, pattern)
    gene_list_neg = gen_uniprot_id_list_neg(gene_list_pos, pattern)
    index = create_index(gene_list_pos, mg)
    fill_gene(gene_list_pos, index, path)
    fill_gene(gene_list_neg, index, path)
    return {"positif": gene_list_pos, "negatif": gene_list_neg}

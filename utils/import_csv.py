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
import pandas as pd
import csv
import os
from biothings_client import get_client
import mygene
import requests
from utils.tools import find_pattern, is_metazoan, print_trace


def import_csv(csv):
    df = pd.read_csv(csv)
    # Convert data into category
    for cat in df.columns:
        if cat != "position":
            df[cat] = df[cat].astype('category')
    return df


def request_gene_id(geneID, s):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1' \
              '&singlecopy=1&limit=1' % geneID
    response = s.get(request)
    if response.status_code == 200:
        return response.json()
    else:
        print("status code for %s = %s" % (request, response.status_code))
        return None


def uniprotid_to_geneid(uniprotid_list):
    mg = mygene.MyGeneInfo()
    if len(uniprotid_list):
        return mg.querymany(uniprotid_list, scope='symbol,accession',
                            fields='uniprot, taxid', species="all", as_dataframe=True)
    else:
        return []


def request_cluster_id(clusterID, path, s):
    name = "%s.fasta" % clusterID
    path2fastas = "%s/fastas" % path
    path2file = "%s/%s" % (path2fastas, name)
    if not os.path.exists(path2fastas):
        os.mkdir(path2fastas)
    if not os.path.exists(path2file):
        request_odb = 'http://www.orthodb.org/fasta?id=%s' % clusterID
        resp = s.get(request_odb)
        if resp.status_code == 200:
            request_api = "curl %s -o %s" % (request_odb, path2file)
            os.system(request_api)
        else:
            print("status code for %s = %s" % (request_odb, resp.status_code))



def import_ortholog(csv_file, pattern):

    print("Parsing csv")

    path = os.path.dirname(os.path.dirname(csv_file))
    file_name = os.path.basename(csv_file)
    index_file = '%s/csv/%s/index_%s_%s.csv' % (path, pattern, file_name[:-4], pattern)

    mt = get_client("taxon")

    print("Extracting %s phosphorylation site" %pattern)
    df = import_csv(csv_file)
    sub_df = df[df["code"] == pattern]
    uniprot_id_list = []

    if os.path.exists(index_file) and os.path.getsize(index_file) > 0:
        index_df = pd.read_csv(index_file, sep=';')
        uniprot_id_list = index_df["uniprotID"].value_counts().keys().tolist()

    print("Preparing queries")
    uniprot_to_convert = set(sub_df["acc"].tolist()) - set(uniprot_id_list)
    resp = uniprotid_to_geneid(uniprot_to_convert)

    with open(index_file, 'a+', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        g.seek(0)
        first_char = g.read(1)
        if not first_char:
            writer.writerow(['uniprotID', 'geneID', 'taxID', 'metazoan', 'code', 'sequence',
                            'pos_sites', "neg_sites", 'clusterID'])

        with requests.Session() as s:
            group_acc_seq = sub_df.groupby(["acc", "sequence"])
            length = len(group_acc_seq)
            for i, ((acc, seq), group)in enumerate(group_acc_seq):
                print_trace(i, length, "import %s from csv file" %acc)
                if acc not in uniprot_id_list:

                    # Info from phospho.ELM
                    pos_list = []
                    for position in group["position"]:
                        if position - 1 not in pos_list:
                            pos_list.append(position - 1)
                    neg_list = []
                    for m in find_pattern(pattern, seq):
                        new_position = round((m.end() + m.start() - 1) / 2)
                        neg = True
                        for pos in pos_list:
                            if abs(new_position - pos) <= 50:
                                neg = False
                        if neg:
                            neg_list.append(new_position)

                    # Info from gene
                    r = resp.loc[[acc]]
                    geneID = None
                    taxID = None

                    if "_id" in r:
                        geneID = r["_id"].values[0]
                    if "taxid" in r:
                        taxID = r["taxid"].values[0]
                    metazoan = is_metazoan(taxID, mt)

                    # Info from orthoDB
                    request = request_gene_id(geneID, s)
                    clusterID = "nan"
                    if request is not None:
                        if len(request["data"]):
                            clusterID = request["data"][0]
                            request_cluster_id(clusterID, path, s)

                    writer.writerow([acc, geneID, taxID, metazoan, pattern, seq, pos_list, neg_list,
                                     clusterID])
    return index_file

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
import ast
import math
import numpy as np
from threading import Thread, RLock
import queue
from source.utils.tools import find_pattern, function_in_thread
from source.utils.sequence import is_metazoan, find_sequence


lock = RLock()


class fill_csv(Thread):
    def __init__(self, group_acc_seq,uniprot_id_list,
                 pattern, mt, path, writer, resp):
        Thread.__init__(self)
        self.group_acc_seq = group_acc_seq
        self.uniprot_id_list = uniprot_id_list
        self.pattern = pattern
        self.mt = mt
        self.path = path
        self.writer = writer
        self.resp = resp

    def run(self):
        for i, (acc, group) in enumerate(self.group_acc_seq):
            print("import %s from csv file" % acc)
            with requests.Session() as s:
                if acc not in self.uniprot_id_list:
                    pos_list = []
                    seq_list = []
                    que = queue.Queue()
                    # Info from gene
                    r = (self.resp).loc[[acc]]
                    geneID = None
                    taxID = None

                    if "_id" in r:
                        geneID = r["_id"].values[0]
                    if "taxid" in r:
                        taxID = r["taxid"].values[0]
                    metazoan = is_metazoan(taxID, self.mt)
                    clusterID = function_in_thread(que, [acc,
                                                         geneID,
                                                         self.path,
                                                         s],
                                                   request_cluster_id)
                    name = "%s.fasta" % acc
                    path2fastas = "%s/fastas" % self.path
                    path2cluster = "%s/%s" % (path2fastas, name)
                    sequence = function_in_thread(que, [path2cluster,
                                                        group["seq_in_window"].tolist(),
                                                        taxID], find_sequence)
                    if sequence != "None":
                        for position, seq in zip(group["position"], group["seq_in_window"]):
                            match = function_in_thread(que, [seq, sequence], find_pattern)
                            tmp_position = None
                            if len(match):
                                for r in match:
                                    if tmp_position is None:
                                        dist = len(sequence) + 1
                                    else:
                                        dist = abs(tmp_position - position)
                                    if dist > abs(6 + r.start() - position):
                                        tmp_position = 6 + r.start()
                                position = tmp_position
                            else:
                                position = None
                            if position not in pos_list and position is not None:
                                pos_list.append(position)
                                seq_list.append(seq)

                        with lock:
                            if len(pos_list):
                                (self.writer).writerow([acc, geneID, taxID, metazoan, self.pattern, seq_list,
                                                        pos_list, clusterID, sequence])


def import_csv(csv):
    df = pd.read_csv(csv, sep="\t", header=None)
    df.columns = [
        'prot_name',
        'acc',
        'position',
        'type',
        'pmids',
        'database',
        'code',
        'tpm',
        'seq_in_window'
    ]
    # Convert data into category
    for cat in df.columns:
        if cat != "position":
            df[cat] = df[cat].astype('category')
    return df


def request_gene_id(geneID, s):
    request = 'http://www.orthodb.org/search?query=%s&ncbi=1' \
              '&singlecopy=1&limit=1&universal=1' % geneID
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


def request_ortho_db(s, id, ncbi, path2file):
    request_odb = 'http://www.orthodb.org/fasta?query=%s&ncbi=%s' % (id, ncbi)
    resp = s.get(request_odb)
    if resp.status_code == 200:
        content = resp.content.decode("ascii")
        try:
            ast.literal_eval(content)
        except:
            request_api = 'wget \'%s\' -O %s' % (request_odb, path2file)
            os.system(request_api)
            return True
    else:
        print("status code for %s = %s" % (request_odb, resp.status_code))
    return False


def request_cluster_id(acc, geneID, path, s):
    cluster_id = "nan"
    name = "%s.fasta" % acc
    path2fastas = "%s/fastas" % path
    path2file = "%s/%s" % (path2fastas, name)
    os.makedirs(path2fastas, exist_ok=True)
    downloaded = False
    if not os.path.exists(path2file):
        downloaded = request_ortho_db(s, acc, 0, path2file)
        if geneID is not None:
            if not downloaded and not math.isnan(float(geneID)):
                downloaded = request_ortho_db(s, geneID, 1, path2file)
    if downloaded or os.path.exists(path2file):
        cluster_id = acc
    return cluster_id


def import_ortholog(csv_file, pattern, nthread):
    print("Parsing csv")

    if os.path.exists("%s/data" % os.path.dirname(csv_file)):
        path = "%s/data" % os.path.dirname(csv_file)
    else:
        path = os.path.dirname(os.path.dirname(csv_file))
    file_name = os.path.basename(csv_file)
    os.makedirs("%s/csv/%s" % (path, pattern), exist_ok=True)
    index_file = '%s/csv/%s/index_%s_%s.csv' % (path, pattern, file_name[:-4], pattern)
    df = import_csv(csv_file)
    mt = get_client("taxon")

    print("Extracting %s phosphorylation site" % pattern)

    uniprot_id_list = []
    if os.path.exists(index_file) and os.path.getsize(index_file) > 0:
        index_df = pd.read_csv(index_file, sep=';')
        uniprot_id_list = index_df["uniprotID"].value_counts().keys().tolist()

    print("Preparing queries")
    uniprot_to_convert = set(df["acc"].tolist()) - set(uniprot_id_list)
    resp = uniprotid_to_geneid(uniprot_to_convert)

    sub_df = df[df["acc"].isin(list(uniprot_to_convert))]
    sub_df.reset_index()
    with open(index_file, 'a+', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        g.seek(0)
        first_char = g.read(1)
        if not first_char:
            writer.writerow(['uniprotID', 'geneID', 'taxID', 'metazoan', 'code',
                             'seq_in_window', 'pos_sites', 'clusterID', 'sequence'])

        group_acc_seq = sub_df.groupby(["acc"], observed=True)
        data_thread = np.array_split(group_acc_seq, nthread)
        thread_list = []
        for data in data_thread:
            thread_list.append(fill_csv(data, uniprot_id_list,
                                        pattern, mt, path, writer, resp))
        for thread in thread_list:
            thread.start()
        for thread in thread_list:
            thread.join()
    return index_file

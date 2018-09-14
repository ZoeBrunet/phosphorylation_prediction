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
import os
import csv
import mygene
import pandas as pd
from xmldict import xml_to_dict


def extract_xml_index(file, prot_set, initial_id_list):
    with open(file) as fd:
        tmp_set = set()
        doc = xml_to_dict(fd.read())
        prot_list = doc["musite"]["protein-list"]["protein"]
        for p in prot_list:
           tmp_set.add(p["accession"])
        return prot_set.union(tmp_set) - initial_id_list


def get_used_prot(used_protein, convert=None):
    initial_id_list = set()
    local_list = []
    index = "%s_index.csv" % os.path.splitext(used_protein)[0]
    if convert is not None:
        convert = pd.read_csv(convert, sep='\t')
        convert.columns = ['uniprotID', 'local_id']
        local_list = convert["local_id"].tolist()
    if os.path.exists(index) and os.path.getsize(index) > 0:
        index_df = pd.read_csv(index, sep=';')
        initial_id_list = set(index_df["initial_id"].value_counts().keys().tolist())

    with open(index, 'a+', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        g.seek(0)
        first_char = g.read(1)
        prot_id = set()
        if not first_char:
            writer.writerow(['uniprotID', 'initial_id'])
        with open(used_protein) as f:
            file_list = f.readlines()
        file_list = [x.strip() for x in file_list]
        for file in file_list:
            extension = os.path.splitext(file)[1]
            if extension == ".txt":
                with open(file) as h:
                    prot = h.readlines()
                    prot_id = prot_id.union(set([x.strip() for x in prot])) - initial_id_list
            if extension == ".xml":
                prot_id = extract_xml_index(file, prot_id, initial_id_list)
        mg = mygene.MyGeneInfo()
        resp = mg.querymany(prot_id, scope='symbol,accession',
                            fields='uniprot', species="all")
        for row in resp:
            prot = row["query"]
            if 'uniprot' in row:
                if "Swiss-Prot" in row["uniprot"]:
                    translation = row["uniprot"]["Swiss-Prot"]
                else:
                    translation = row["uniprot"]["TrEMBL"]
                if type(translation) is not list:
                    translation = [translation]
                for t in translation:
                    writer.writerow([t, prot])
            else:
                if convert is not None:
                    tmp = prot.split(".")[0]
                    if tmp in local_list:
                        writer.writerow([convert[convert["local_id"] == tmp]["uniprotID"].values[0], prot])
                    else:
                        writer.writerow([None, prot])
                else:
                    writer.writerow([None, prot])
            print("add %s" % prot)
    return index


def split_dataset(used_protein, convert, dataset_file):
    benchmark_file = "%s_benchmark.csv" % os.path.splitext(dataset_file)[0]
    training_file = "%s_training.csv" % os.path.splitext(dataset_file)[0]
    index = pd.read_csv(get_used_prot(used_protein, convert), sep=';')
    initial_id_list = set(index["initial_id"].value_counts().keys().tolist())
    uniprotID_list = set(index["uniprotID"].value_counts().keys().tolist())
    dataset = pd.read_csv(dataset_file, sep=';')
    training = dataset[(dataset["uniprotID"].isin(uniprotID_list) |
                        dataset["uniprotID"].isin(initial_id_list) |
                        dataset["geneID"].isin(uniprotID_list) |
                        dataset["geneID"].isin(initial_id_list) |
                        dataset["sequence"].isnull())]
    benchmark = dataset[~dataset.index.isin(training.index)]
    training.to_csv(training_file, index=False, sep=';')
    benchmark.to_csv(benchmark_file, index=False, sep=';')
    return {"benchmark": benchmark_file, "training": training_file}

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
import pandas as pd
import csv


def get_info(file):
    path = os.path.dirname(file)
    file_name = "info_%s" % os.path.basename(file)
    output_file = "%s/%s" % (path, file_name)
    df = pd.read_csv(file, sep=";")
    if not os.path.exists(output_file):
        os.mknod(output_file)
    with open(output_file, 'w', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        writer.writerow(["species", "nb_prot", "nb_sites", "nb_sites_pos", "nb_sites_neg"])
        taxID_list = df["taxID"].value_counts().keys().tolist()
        for taxID in taxID_list:
            sub_df = df[df["taxID"] == taxID]
            nb_prot = set(sub_df["uniprotID"].tolist()).__len__()
            nb_pos = len(sub_df[sub_df["phosphorylation_site"]])
            nb_neg = len(sub_df[~sub_df["phosphorylation_site"]])
            writer.writerow([taxID, nb_prot, len(sub_df), nb_pos, nb_neg])

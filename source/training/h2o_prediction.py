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

import h2o
from h2o.utils.shared_utils import _locate
import csv
import os
import pandas as pd
from pathlib import Path


def h2o_prediction(input_file, prediction, benchmark):
    if not os.path.exists(prediction):
        model_name = os.path.basename(input_file)
        h2o.init(nthreads=-1)
        model = h2o.load_model(input_file)

        df = h2o.import_file(path=_locate(benchmark))

        for cat in ["taxID", "phosphorylation_site", "metazoa"]:
            df[cat] = df[cat].asfactor()
        col_with_nan = []
        for column in df.columns:
            if "freq" in column or "ACH" in column or "IC" in column or "shanon_entropy" in column:
                col_with_nan.append(column)
        for cat in [["geneID", "taxID", "clusterID", "metazoa"] + col_with_nan]:
            df[df[cat] == "nan", cat] = float('nan')

        # Variable selection

        df_names_x = df.names[:]
        for col in ["phosphorylation_site", "uniprotID", "geneID", "sequence", "position", "seq_in_window",
                    "clusterID", "nb_orthologs", "taxID", "nb_orthologs_metazoa", "nb_orthologs_non_metazoa"]:
            df_names_x.remove(col)
        predict = model.predict(df).as_data_frame(use_pandas=True, header=True)
        bench = pd.read_csv(benchmark, sep=';')
        Path(prediction).touch(exist_ok=True)
        with open(prediction, 'r+', newline='') as g:
            writer = csv.writer(g, delimiter=";")
            writer.writerow(["uniprotID", "position", "phosphorylation_site",
                             "sequence", "prediction_%s" % model_name])
            for score, uniprotID, position, \
                phosphorylation_site, seq in zip(predict["True"],
                                                 bench["uniprotID"],
                                                 bench["position"],
                                                 bench["phosphorylation_site"],
                                                 bench["sequence"]):
                writer.writerow([uniprotID, position, phosphorylation_site, seq, score])

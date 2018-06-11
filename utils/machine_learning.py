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
import sys
import h2o
import os
from h2o.utils.shared_utils import _locate
from h2o.automl import H2OAutoML


def train_model(file, max_mod):
    h2o.init(nthreads=-1, max_mem_size="1g")

    print("Import and Parse data")
    df = h2o.import_file(path=_locate(file))

    # Initialisation

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
    for col in ["phosphorylation_site", "uniprotID", "geneID",
                "clusterID", "nb_orthologs"]:
         df_names_x.remove(col)

    # Model creation

    aml = H2OAutoML(max_models=max_mod, seed=1)
    aml.train(x=df_names_x, y="phosphorylation_site",
              training_frame=df)
    lb = aml.leaderboard
    path = os.path.dirname(os.path.dirname(os.path.dirname(file)))
    pattern = os.path.dirname(file)[-1:]

    # Save models
    path2models = "%s/models/%s/%smodels" % (path, pattern, max_mod)
    for id in list(lb['model_id'].as_data_frame().iloc[:, 0]):
        model = h2o.get_model(id)
        h2o.save_model(model, path=path2models,
                       force=True)
    with open('%s/info.txt' % path2models, 'w', newline='') as g:
        orig_stdout = sys.stdout
        sys.stdout = g
        print(lb)
        sys.stdout = orig_stdout

    h2o.cluster().shutdown

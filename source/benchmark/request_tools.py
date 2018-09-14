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
import mechanicalsoup
import requests
import io
import csv
import os
import math
import numpy as np
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
from functools import reduce
from pathlib import Path
from requests.exceptions import ConnectionError
from matplotlib.pyplot import cm
from source.training.h2o_prediction import h2o_prediction


def get_statistics(TP, TN, FP, FN):
    recall = TP/(TP + FN) if (TP + FN) != 0 else None
    precision = TP / (TP + FP) if (TP + FP) != 0 else None
    specificity = TN / (TN + FP) if (TN + FP) != 0 else None
    antispecificity = 1 - specificity if specificity is not None else None
    MCC = ((TP * TN)-(FP * FN))/math.sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)) \
        if (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN) != 0 else None
    F1 = 2*(precision * recall)/(precision + recall) if ((precision is not None and recall is not None) and (int(precision + recall) != 0)) else None
    return {"recall": recall, "precision": precision,
            "specificity": specificity, "antispecificity": antispecificity,
            "MCC": MCC, "F1": F1}


def get_confusion_matrix(input_file, threshold, label):
    df = pd.read_csv(input_file, sep=";")
    df
    TP = TN = FP = FN = 0
    for predict, phospho in zip(df[label], df["phosphorylation_site"]):
        if predict <= threshold and phospho == True:
            FN += 1
        if predict <= threshold and phospho == False:
            TN += 1
        if predict > threshold and phospho == True:
            TP += 1
        if predict > threshold and phospho == False:
            FP += 1
    if len(df) != int(TP + TN + FP + FN):
        print("Problem occured during confusion matrix generation")
    return{"TP": TP, "TN": TN, "FP": FP, "FN": FN}


def RF_Phos(fasta, pos):
    browser = mechanicalsoup.StatefulBrowser()
    browser.open("http://bcb.ncat.edu/RF_Phos/")
    browser.select_form('form[id="SeqForm"]')
    browser["fastaSeq"] = fasta
    browser["serine"] = "serine"
    browser["threonine"] = "threonine"
    browser["tyrosine"] = "tyrosine"
    response = browser.submit_selected()
    html = response.text
    soup = BeautifulSoup(html)
    table_S = soup.find("table", attrs={"id": "tb_details_id_S"})
    table_T = soup.find("table", attrs={"id": "tb_details_id_T"})
    table_Y = soup.find("table", attrs={"id": "tb_details_id_Y"})
    new_header = None
    df = None
    if table_S is not None:
        new_header = pd.read_html(str(table_S))[0].iloc[0]
    if table_T is not None:
        new_header = pd.read_html(str(table_T))[0].iloc[0]
    if table_Y is not None:
        new_header = pd.read_html(str(table_Y))[0].iloc[0]
    if table_S is not None and table_T is not None and table_Y is not None:
        frame = [pd.read_html(str(table_S))[0][1:], pd.read_html(str(table_T))[0][1:],
                 pd.read_html(str(table_Y))[0][1:]]
        df = pd.concat(frame)
    if new_header is not None:
        df.columns = new_header
    if df is not None:
        for position, score in zip(df["Position"], df["Score"]):
            if int(position) == pos:
                return score
    return None


def musite(fasta, pos):
    seq = fasta.split("\n")[1]
    request_S_T = "http://musite.net/service?type=pred&seq=%s&model=Phosphorylation.Eukaryote.ser.thr&sp=.0" % seq
    request_Y = "http://musite.net/service?type=pred&seq=%s&model=Phosphorylation.Eukaryote.tyr&sp=.0" % seq
    try:
        resp_S_T = requests.get(request_S_T)
    except ConnectionError:
        resp_S_T = requests.get(request_S_T)
    try:
        resp_Y = requests.get(request_Y)
    except ConnectionError:
        resp_Y = requests.get(request_Y)
    df_S_T = pd.read_csv(io.StringIO(resp_S_T.content.decode("utf-8")),
                         sep="\t", names=["Residue", "Position", "Score"])
    df_Y = pd.read_csv(io.StringIO(resp_Y.content.decode("utf-8")),
                       sep="\t", names=["Residue", "Position", "Score"])
    df = df_S_T if pos in df_S_T["Position"].tolist() else df_Y
    if len(df[df["Position"] == pos]):
        return df[df["Position"] == pos]["Score"].values[0]


def tools_prediction(inputfile, tool_list):
    path = "%s/prediction" % os.path.dirname(inputfile)
    os.makedirs(path, exist_ok=True)
    benchmark = pd.read_csv(inputfile, sep=';')
    gb = benchmark.groupby(["uniprotID", "position", "phosphorylation_site", "sequence"]).groups
    for key, value in tool_list.items():
        if value:
            outputfile = "%s/%s.csv" % (path, key)
            formergb = []
            if os.path.exists(outputfile) and os.path.getsize(outputfile) > 0:
                formerpos = pd.read_csv(outputfile, sep=';')
                formergb = formerpos.groupby(("uniprotID", "position",
                                              "phosphorylation_site", "sequence")).groups
            gbtoconvert = list(set(gb) - set(formergb))
            l = len(gbtoconvert)
            with open(outputfile, 'a+', newline='') as g:
                writer = csv.writer(g, delimiter=";")
                g.seek(0)
                first_char = g.read(1)
                if not first_char:
                    writer.writerow(["uniprotID", "position", "phosphorylation_site",
                                     "sequence", "prediction_%s" % key])
                for i, g in enumerate(gbtoconvert):
                    uniprotID, position, phosphorylation_site, seq = g
                    pos = position + 1
                    fasta = ">\n%s" % seq
                    prediction = eval("%s(fasta, pos)" % key)
                    writer.writerow([uniprotID, position, phosphorylation_site, seq, prediction])
                    print("importation for %s %s/%s = %s" % (key,
                                                             str(i + 1), str(l),
                                                             str(round(((i + 1) / l) * 100, 2)) + "%"))
    return path


def model_prediction(output_directory, model_directory, models, benchmark):
    if not len(models):
        model_family = ["StackedEnsemble_AllModels",
                        "StackedEnsemble_BestOfFamily",
                        "GBM",
                        "XRT",
                        "DRF",
                        "DeepLearning",
                        "GLM"]
        info = "%s/info.txt" % model_directory
        with open(info) as g:
            txt = g.readlines()
            lines = [x.strip() for x in txt]
        for line in lines:
            model_name = line.split(" ")[0]
            for mod in model_family:
                if mod in model_name:
                    models.append(model_name)
                    model_family.remove(mod)
    for model in models:
        output_file = "%s/%s.csv" % (output_directory, model)
        h2o_prediction("%s/%s" % (model_directory, model), output_file, benchmark)
    return models


def benchmark_result(path, tool_list):
    prediction_df = []
    column_prediction = []
    result = "%s/global_results.csv" % path
    for key in tool_list:
        path2pred = "%s/%s.csv" % (path, key)
        prediction_df.append(pd.read_csv(path2pred, sep=';'))
        column_prediction.append("prediction_%s" % key)
    global_df = reduce(lambda left, right: pd.merge(left, right, on=['uniprotID', "position",
                                                                     "phosphorylation_site", "sequence"]),
                       prediction_df)
    global_df = global_df.dropna()
    df_pos = global_df[global_df["phosphorylation_site"]]
    df_neg = global_df[~global_df["phosphorylation_site"]].sample(n=len(df_pos))
    frames = [df_pos, df_neg]
    merged = pd.concat(frames)
    merged.to_csv(result, index=False, sep=';')
    return result


def get_result_info(result_file, step):
    stat_file = "%s/stat.csv" % os.path.dirname(result_file)
    best_stat_file = "%s/best_stat.csv" % os.path.dirname(result_file)
    thresholds = [x / step for x in range(1, step)]
    result = pd.read_csv(result_file, sep=';')
    pred = []
    for column in result.columns:
        if "prediction" in column:
            pred.append(column)
    Path(stat_file).touch(exist_ok=True)
    with open(stat_file, 'r+', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        writer.writerow(["tool", "threshold", "TP", "TN", "FP", "FN",
                         "recall", "precision", "specificity",
                         "antispecificity", "MCC", "F1"])
        for col in pred:
            for threshold in thresholds:
                confusion_matrix = get_confusion_matrix(result_file, threshold, col)
                TP = confusion_matrix["TP"]
                TN = confusion_matrix["TN"]
                FP = confusion_matrix["FP"]
                FN = confusion_matrix["FN"]
                stat = get_statistics(TP, TN, FP, FN)
                recall = stat["recall"]
                precision = stat["precision"]
                specificity = stat["specificity"]
                antispecificity = stat["antispecificity"]
                MCC = stat["MCC"]
                F1 = stat["F1"]
                writer.writerow([col, threshold, TP, TN, FP, FN,
                                 recall, precision, specificity,
                                 antispecificity, MCC, F1])
    frames = []
    stat = pd.read_csv(stat_file, sep=';')
    for col in pred:
        tmp = stat[stat["tool"] == col].dropna()
        max = tmp['MCC'].max()
        frames.append(tmp[tmp['MCC'] == max].head(1))
    best_stat = pd.concat(frames)
    best_stat.to_csv(best_stat_file, index=False, sep=';')
    return stat_file


def plot_curves(stat_file):
    stat = pd.read_csv(stat_file, sep=';')
    tools = stat["tool"].value_counts().keys().tolist()
    plt.subplot(121)
    plt.xlabel('recall')
    plt.ylabel('precision')
    plt.title('Precision-Recall Curve')
    color = iter(cm.rainbow(np.linspace(0, 1, len(tools))))
    for tool in tools:
        c = next(color)
        recall = stat[stat["tool"] == tool].sort_values(by=['threshold'])["recall"]
        precision = stat[stat["tool"] == tool].sort_values(by=['threshold'])["precision"]
        plt.plot(recall, precision, c=c, label=tool.replace('prediction_', ''))
    plt.legend(loc="lower left")
    plt.subplot(122)
    plt.xlabel('1 - specificity')
    plt.ylabel('recall')
    plt.title('ROC Curve')
    colorbis = iter(cm.rainbow(np.linspace(0, 1, len(tools))))
    for tool in tools:
        c = next(colorbis)
        recall = stat[stat["tool"] == tool].sort_values(by=['threshold'])["recall"]
        antispecificity = stat[stat["tool"] == tool].sort_values(by=['threshold'])["antispecificity"]
        antispecificity_without_none = [0 if elem is None else elem for elem in antispecificity]
        auc = abs(np.trapz(recall, antispecificity_without_none))
        plt.plot(antispecificity_without_none, recall, c=c, label='%s (area = %0.4f)' %
                                                                  (tool.replace('prediction_', ''), auc))
    lims = [np.min([0, 0]), np.max([1, 1]), ]
    plt.plot(lims, lims, '--b', label='Random classifier')
    plt.legend(loc="lower right")
    plt.show()

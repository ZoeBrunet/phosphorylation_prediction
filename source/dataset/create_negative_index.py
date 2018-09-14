import pandas as pd
import csv
import ast
import re
import os
from source.utils.window import create_window


def write_neg_sites(output_file, input_file, file_list, pattern):
    if not os.path.exists(output_file):
        genes = pd.read_csv(input_file, sep=";")
        df_list = []
        sequence_list = []
        for file in file_list:
            df = pd.read_csv(file, sep=";")
            df_list.append(df)
            sequence_list.append(df["sequence"].tolist())
        length = len(genes)

        with open(output_file, 'a+', newline='') as g:
            writer = csv.writer(g, delimiter=";")
            writer.writerow(['uniprotID', 'geneID', 'taxID', 'metazoan', 'code',
                             'neg_seq_in_window', 'neg_sites', 'clusterID', 'sequence'])
            for i, (uniprotID, geneID, taxID, metazoan, code, \
                    pos_sites, clusterID, \
                    sequence) in enumerate(zip(genes['uniprotID'], genes['geneID'], genes['taxID'], genes['metazoan'], \
                                               genes['code'], genes['pos_sites'], \
                                               genes['clusterID'], genes['sequence'])):
                position_list = ast.literal_eval(pos_sites)
                for df, seq in zip(df_list, sequence_list):
                    if sequence in seq and str(sequence) != "nan":
                        position_list += ast.literal_eval(df[df["sequence"] == sequence]["pos_sites"].values[0])
                p = re.compile(pattern)
                result = []
                neg_sites = []
                seq_in_window =[]
                if str(sequence) != "nan":
                    for m in p.finditer(sequence):
                        result.append(m)
                    for m in result:
                        new_position = round((m.end() + m.start() - 1) / 2)
                        neg = True
                        for position in position_list:
                            if abs(new_position - position) <= 50:
                                neg = False
                        if neg:
                            neg_sites.append(new_position)
                            rel_window = create_window(sequence, new_position, 13, True)
                            seq_in_window.append(sequence[rel_window[0][0]:rel_window[1][1] + 1])
                print("Import negative site %s/%s = %s" % (str(i + 1), str(length),
                                                           str(round(((i + 1) / length) * 100, 2)) + "%"))
                writer.writerow([uniprotID, geneID, taxID, metazoan, code,
                                 seq_in_window, neg_sites, clusterID, sequence])
    return output_file


def merge_index(file1, file2, output):
    if not os.path.exists(output):
        a = pd.read_csv(file1, sep=";")
        b = pd.read_csv(file2, sep=";")
        merged = a.merge(b, on=['uniprotID', 'geneID', 'taxID', 'metazoan', 'code', 'clusterID', 'sequence'])
        merged.to_csv(output, index=False, sep=';')
    return output


def write_final_index(input_file, file_list, pattern):
    name = "index_neg_%s.csv" % pattern
    path = os.path.dirname(input_file)
    output_file = "%s/%s" % (path, name)
    output_file = write_neg_sites(output_file, input_file, file_list, pattern)
    output = "%s/final_index_%s.csv" % (path, pattern)
    return merge_index(input_file, output_file, output)


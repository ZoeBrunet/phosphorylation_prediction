import os
import pandas as pd


def filter_species(file, taxID_list):
    path = os.path.dirname(file)
    suffix = ""
    df = pd.read_csv(file, sep=";")
    for t in taxID_list:
        suffix += str("_%s" % t)
    file_name = "%s%s.csv" % (os.path.basename(file)[:-4], suffix)
    path2file = "%s/%s" % (path, file_name)
    if not os.path.exists(path2file):
        df = df[df["taxID"].isin(taxID_list)]
        df.to_csv(path2file, index=False, sep=';')
    return path2file

import os
import pandas as pd


def filter_species(file, taxID_list, suffix):
    df = pd.read_csv(file, sep=";")
    path2file = "%s%s.csv" % (os.path.splitext(file)[0], suffix)
    if not os.path.exists(path2file):
        df = df[df["taxID"].isin(taxID_list)]
        df.to_csv(path2file, index=False, sep=';')
    return path2file

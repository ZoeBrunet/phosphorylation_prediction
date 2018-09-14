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
import sys
import pandas as pd
import h2o
from pprint import pprint
import matplotlib.pyplot as plt


def plot_variable_distribution(color_list, ax, dict, str_contain, title):
    length = 0
    for (key, value), color_curve in zip(dict.items(), color_list):
        plt.xlabel('feature')
        plt.ylabel("%")
        plt.title(title)
        former_length = length
        value = value[value["Variable"].str.contains(str_contain)].sort_values(by=['Variable'])
        length += len(value)
        percentage = value["Percentage"].tolist()
        feature = [x for x in range(former_length, length)]
        axe = [-0.003] * (length - former_length)
        plt.plot(feature, percentage, '-%s' % color_curve)
        plt.plot(feature, axe, '-%s' % color_curve, linewidth=3)
        ax.annotate(key, xy=(former_length, -0.002), xytext=(former_length, -0.002),  color=color_curve)


def plot_varimp(model, path2info, path2varimp):
    h2o.init(nthreads=-1)
    mod = h2o.load_model(model)
    var_df = pd.DataFrame(mod.varimp(),
                          columns=["Variable", "Relative Importance", "Scaled Importance", "Percentage"])
    feature_list = ["ACH", "IC", "freq", "entropy"]
    color_list = ['r', 'g', 'c', 'k']
    dict = {}
    for var_string in feature_list:
        dict[var_string] = var_df[var_df["Variable"].str.contains(var_string)].sort_values(by=['Variable'])

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    plot_variable_distribution(color_list, ax1, dict, "", "Variable importance")

    ax2 = fig.add_subplot(223)
    plot_variable_distribution(color_list, ax2, dict, "_metazoa_", "Metazoa variable importance")

    ax3 = fig.add_subplot(224)
    plot_variable_distribution(color_list, ax3, dict, "_nonmetazoa_", "Non_metazoa variable importance")

    var_df.to_csv(path2varimp, index=False, sep=';')
    with open(path2info, 'w', newline='') as g:
        orig_stdout = sys.stdout
        sys.stdout = g
        pprint(mod.get_params())
        sys.stdout = orig_stdout

    plt.show()

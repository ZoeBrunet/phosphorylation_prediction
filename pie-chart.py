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
import plotly
import plotly.graph_objs as go
from source.utils.parser import boxplot_parser

args = boxplot_parser(sys.argv[1:])
files = [str(item) for item in args.files.split(',')]
feature = args.feature
filename = args.filename

for file in files:
    df = pd.read_csv(file, sep=";")
    label = df[feature].value_counts().keys().tolist()
    values = df[feature].value_counts().tolist()
    trace = go.Pie(labels=label, values=values, textinfo='value')
    plotly.offline.plot([trace], filename=filename, image='png')

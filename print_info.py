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
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import sys
from utils.import_csv import import_csv
from utils.parser import info_parser


def print_pie(csv, column, username, apikey, caption):
    df = import_csv(csv)
    plotly.tools.set_credentials_file(username=username, api_key=apikey)
    label = df[column].value_counts().keys().tolist()
    values = df[column].value_counts().tolist()
    trace = go.Pie(labels=label, values=values, textinfo='value')
    py.plot([trace], filename=caption)


args = info_parser(sys.argv[1:])
print_pie(args.file, args.column, args.username, args.apikey, args.caption)

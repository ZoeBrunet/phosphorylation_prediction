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
plotly.tools.set_credentials_file(username='usernam', api_key='apikey')

import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go


def import_csv(csv):
    df = pd.read_csv(csv)
    # Convert data into category
    for cat in df.columns:
        if cat != "position":
            df[cat] = df[cat].astype('category')
    return df


def print_pie(df):
    label = df['code'].value_counts().keys().tolist()
    values = df['code'].value_counts().tolist()
    trace = go.Pie(labels=label, values=values, textinfo='value')
    py.plot([trace], filename='Nature_phosphorylation_sites')


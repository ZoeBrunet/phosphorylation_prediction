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
from source.utils.parser import compare_tools_parser
from source.benchmark.request_tools import *

args = compare_tools_parser(sys.argv[1:])
benchmark = args.benchmark
tool_dict = {"musite": args.musite, "RF_Phos": args.rfp}
models = args.models
if models is not None:
    models = [str(item) for item in models.split(',')]
else:
    models = []

directory = tools_prediction(benchmark, tool_dict)
model_list = model_prediction(directory, args.model_directory, models, benchmark)
tool_list = [key for (key, value) in tool_dict.items() if value]
model_list += tool_list
results = benchmark_result(directory, model_list)
stat = get_result_info(results, args.step)
plot_curves(stat, model_list)

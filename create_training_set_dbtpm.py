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
from utils.parser import dataset_parser
from utils.create_dataset import create_training_set
from utils.find_neg_sites import *


args = dataset_parser(sys.argv[1:])
input_file = [str(item) for item in args.file.split(',')]
patterns = [str(item) for item in args.pattern.split(',')]
if len(patterns) != len(input_file):
    patterns = [patterns[0]] * len(input_file)
dic = {}
files = {}
dic["phospho_sites"] = {}

for file, pattern in zip(input_file, patterns):
    output_file = None
    if pattern in files:
        output_file = str(files[str(pattern)])
    dic = create_training_set(pattern, file, args.max_window,
                              args.nthread, dic["phospho_sites"], phospho_ELM=False, color=args.color,
                              align_ortho_window=args.ortholog, output_file=output_file)
    if pattern not in files:
        files[str(pattern)] = dic["file"]
pos_sites = dic["phospho_sites"]
for pattern in list(set(patterns)):
    find_neg(pos_sites, pattern, args.max_window, args.color, args.nthread, files[pattern],
             args.ortholog)

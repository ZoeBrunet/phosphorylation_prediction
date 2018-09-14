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
from source.dataset.import_csv import import_ortholog
from source.dataset.filter_index import remove_redundancy
from source.dataset.create_negative_index import write_final_index
from source.utils.parser import dataset_parser
from source.dataset.create_dataset import create_training_set
from source.utils.filter_species import filter_species

args = dataset_parser(sys.argv[1:], 'Run create_training set to get the scores of datasets')
input_file = [str(item) for item in args.file.split(',')]
patterns = [str(item) for item in args.pattern.split(',')]
if len(patterns) != len(input_file):
    patterns = [patterns[0]] * len(input_file)
species = []
if args.species is not None:
    species = [int(item) for item in args.species.split(',')]

dic = {}
files = {}
file_list = []
index_list = []

for file, pattern in zip(input_file, patterns):
    index_file = import_ortholog(file, pattern, args.nthread)
    index_file = filter_species(index_file, species)
    index_file = remove_redundancy(index_file)
    file_list.append(index_file)

for file, pattern in zip(file_list, patterns):
    tmp = file_list.copy()
    tmp.remove(file)
    index_list.append(write_final_index(file, tmp, pattern))

for file, pattern in zip(index_list, patterns):
    output_file = None
    if pattern in files:
        output_file = str(files[str(pattern)])
    create_training_set(pattern, args.max_window,
                        args.nthread, file, phospho_ELM=False, color=args.color,
                        align_ortho_window=args.ortholog, output_file=output_file)

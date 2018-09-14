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
from source.utils.parser import dataset_parser
from source.dataset.create_dataset import create_training_set

args = dataset_parser(sys.argv[1:], 'Run final_index_to_dataset to get dataset from index')
file = args.file
pattern = args.pattern

create_training_set(pattern, args.max_window,
                    args.nthread, file, phospho_ELM=False, color=args.color,
                    align_ortho_window=args.ortholog)

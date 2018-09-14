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
from source.utils.parser import split_parser
from source.benchmark.list_forbiden_genes import split_dataset
from source.utils.get_info_on_dataset import get_info

args = split_parser(sys.argv[1:])
dict = split_dataset(args.used_protein, args.convert, args.dataset)

for value in dict.items():
    get_info(value[1])
# Print les infos pour chaque dataset

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


import os
import argparse
import sys
from convert_id import import_ortholog
from utils import scoring


def parse_args(args):
    parser = argparse.ArgumentParser(description='Run align-file to create'
                                                 'dataset from phospho.ELM dump')
    parser.add_argument('path',
                        help='Where is your folder ?')
    parser.add_argument('file',
                        help='Input file containing examples')
    return parser.parse_args(args)


def align_file(path, file_name):
    import_ortholog("%s/%s" %(path, file_name))
    string = "S"
    path2fastas = '%s/fastas' % path
    for file in os.listdir(path2fastas):
        max_window = 15
        path2file = '%s/%s' % (path2fastas, file)
        print(scoring(string, path2file, max_window))
        remove_useless_file = "rm %s" % path2file
        os.system(remove_useless_file)


args = parse_args(sys.argv[1:])
align_file(args.path, args.file)

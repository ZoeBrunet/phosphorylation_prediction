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
import argparse


def arg_parser(parser):
    parser.add_argument('file',
                        help='Input absolute path file containing examples')
    parser.add_argument('max_window', type=int, nargs='?', default=15,
                        help='Size of the windows which contain your pattern')


def common_parser(args, desc):
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('pattern',
                        help='Input Python regular expression you want to detect')
    arg_parser(parser)
    return parser.parse_args(args)


def IC_parser(args):
    parser = argparse.ArgumentParser(description='Run information_content '
                                                 'to get information content')
    arg_parser(parser)
    return parser.parse_args(args)
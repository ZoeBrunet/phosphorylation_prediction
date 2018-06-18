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


def file_parser(parser):
    parser.add_argument('file',
                        help='Input absolute path file containing examples')


def arg_parser(parser):
    file_parser(parser)
    parser.add_argument('max_window', type=int, nargs='?', default=15,
                        help='Size of the windows which contain your pattern')


def dataset_parser(args):
    parser = argparse.ArgumentParser(description='Run program to get training set in csv files')
    parser.add_argument('pattern',
                        help='Input Python regular expression you want to detect')
    arg_parser(parser)
    parser.add_argument('--progression', action='store_true',
                        help='Enable this bool to display progression')
    parser.add_argument('--color', action='store_true',
                        help='Enable this bool to have color in output console')
    parser.add_argument('--ortholog', action='store_true',
                        help='Enable this bool to display ortholog in the window')
    return parser.parse_args(args)


def IC_parser(args):
    parser = argparse.ArgumentParser(description='Run information_content '
                                                 'to get information content')
    arg_parser(parser)
    return parser.parse_args(args)


def muscle_parser(args):
    parser = argparse.ArgumentParser(description='Run run_muscle '
                                                 'to align orthologs from fasta file')
    file_parser(parser)
    return parser.parse_args(args)


def shanon_parser(args):
    parser = argparse.ArgumentParser(description='Run shanon_entropy to get the '
                                                 'shanon entropy of each position '
                                                 'in the alignment')
    file_parser(parser)
    return parser.parse_args(args)


def info_parser(args):
    parser = argparse.ArgumentParser(description='Run print_info to plot pie chart ')
    file_parser(parser)
    parser.add_argument('column',
                        help='Input the name of the column which interess you')
    parser.add_argument('caption',
                        help='Input the title of the figure')
    return parser.parse_args(args)


def enrichment_parser(args):
    parser = argparse.ArgumentParser(description='Run enrich_csv '
                                                 'to merge your csv files')
    parser.add_argument('csv1',
                        help='Input absolute path file containing your'
                             'first csv')
    parser.add_argument('csv2',
                        help='Input absolute path file containing your'
                             'second csv')
    return parser.parse_args(args)


def ACH_parser(args):
    parser = argparse.ArgumentParser(description='Run ACH to get the average '
                                                 'cumulative hydrophobicity of a sequence')
    parser.add_argument('seq',
                        help='Input sequence of amino acid')
    return parser.parse_args(args)


def machine_learning_parser(args):
    parser = argparse.ArgumentParser(description='Run machine learning to create '
                                                 'model from dataset')
    file_parser(parser)
    parser.add_argument('-max_models', default=None,
                        help='Input number of model you want to create')
    parser.add_argument('-max_time', default=None,
                        help='Input the time in second to rum autoML')
    parser.add_argument('-max_mem_size', default='1g',
                        help='Input the maximum size, in bytes, of the memory allocation pool to H2O. '
                             'This value must a multiple of 1024 greater than 2MB. '
                             'Append the letter m or M to indicate megabytes, '
                             'or g or G to indicate gigabytes.')
    return parser.parse_args(args)

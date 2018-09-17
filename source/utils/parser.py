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


def dataset_parser(args, desc):
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('pattern',
                        help='Input Python regular expression you want to detect')
    parser.add_argument('--nthread', type=int, nargs='?', default=1,
                        help='Input number of thread you want to use')
    parser.add_argument('--species', nargs='?', default=None,
                        help='Input taxid list coma separated')
    arg_parser(parser)
    parser.add_argument('--color', action='store_true',
                        help='Enable this bool to have color in output console')
    parser.add_argument('--ortholog', action='store_true',
                        help='Enable this bool to display ortholog in the window')
    return parser.parse_args(args)


def score_parser(args, desc):
    parser = argparse.ArgumentParser(description=desc)
    arg_parser(parser)
    return parser.parse_args(args)


def muscle_parser(args):
    parser = argparse.ArgumentParser(description='Run run_muscle '
                                                 'to align orthologs from fasta file')
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
    parser.add_argument('directory_name',
                        help='Input the name of the directory where you want to store models')
    parser.add_argument('-max_models', default=None,
                        help='Input number of model you want to create')
    parser.add_argument('-max_time', default=None,
                        help='Input the time in second to rum autoML')
    parser.add_argument('-max_mem_size', default='2g',
                        help='Input the maximum size, in bytes, of the memory allocation pool to H2O. '
                             'This value must a multiple of 1024 greater than 2MB. '
                             'Append the letter m or M to indicate megabytes, '
                             'or g or G to indicate gigabytes.')
    return parser.parse_args(args)


def boxplot_parser(args):
    parser = argparse.ArgumentParser(description='Run program to get boxplot or pie-chart')
    parser.add_argument('feature',
                        help='Input the feature you want to plot')
    parser.add_argument('files',
                        help='Input list of csv with "," as separator')
    parser.add_argument('filename',
                        help='Input the path of the output')
    parser.add_argument('-names', default=None,
                        help='Optional input list of name "," as separator for boxplot')
    return parser.parse_args(args)


def model_parser(args):
    parser = argparse.ArgumentParser(description='Run program to get info on model')
    parser.add_argument('model', help='Input path to model')
    return parser.parse_args(args)


def split_parser(args):
    parser = argparse.ArgumentParser(description='Run program to split dataset '
                                                 'into validation and training set')
    parser.add_argument('dataset', help='Input path to dataset')
    parser.add_argument('used_protein', help='Input path to file in which are '
                                             'Musite and RF-phos dataset')
    parser.add_argument('-convert', default=None,
                        help='Optional input convert file to translate id into uniprotid')
    return parser.parse_args(args)


def compare_tools_parser(args):
    parser = argparse.ArgumentParser(description='Run program to compare models predictions')
    parser.add_argument('benchmark', help='Input path to benchmark')
    parser.add_argument('model_directory', help='Input the directory in whom models are stored')
    parser.add_argument('--musite', action='store_true',
                        help='Enable this bool to enable Musite prediction')
    parser.add_argument('--rfp', action='store_true',
                        help='Enable this bool to enable RF-Phos prediction')
    parser.add_argument('-step', default=100,
                        help='Optional input step for threshold')
    parser.add_argument('-models', default=None,
                        help='Optional input list of models you want to test')
    return parser.parse_args(args)

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
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC


def get_align_info(file):
    alignment = AlignIO.read(open(file), "fasta", alphabet=IUPAC.extended_protein)
    return AlignInfo.SummaryInfo(alignment)


def get_pssm(summary_align):
    consensus = summary_align.dumb_consensus()
    my_pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore=['-'])
    return my_pssm


def get_information_content(file, window):
    summary_align = get_align_info(file)
    get_pssm(summary_align)
    info_content = summary_align.information_content(window[0], window[1],
                                                     log_base=10,
                                                     chars_to_ignore=['-'])
    return info_content


def parse_args_ic(args):
    parser = argparse.ArgumentParser(description='Run scoring to show '
                                                 'the information content')
    parser.add_argument('file',
                        help='Input file containing examples')
    parser.add_argument('pos1', type=int,
                        help='Begin of the window')
    parser.add_argument('pos2', type=int,
                        help='End of the window')
    return parser.parse_args(args)


args = parse_args_ic(sys.argv[1:])
window = [args.pos1, args.pos2]
print(get_information_content(args.file, window))

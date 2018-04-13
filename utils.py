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
import re
import argparse
from io import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


def parse_args(args):
    parser = argparse.ArgumentParser(description='Run scoring to detect '
                                                 'the frequency of some pattern')
    parser.add_argument('pattern',
                        help='Input Python regular expression you want to detect')
    parser.add_argument('file',
                        help='Input file containing examples')
    parser.add_argument('max_window', type=int, nargs='?', default=15,
                        help='Size of the windows which contain your pattern')
    return parser.parse_args(args)


def find_pattern(pattern, seq):
    p = re.compile(pattern)
    result = []
    for m in p.finditer(seq):
        result.append(m)
    return result


def fill_score_table(score, m, align):
    window_length = m.end() - m.start()
    upper = round((window_length / 2))
    if window_length > 0:
        for i in range(0, upper + 1):
            coef = (i + 1) / (align.__len__() * window_length)
            score[m.start() + i] += coef
            if m.end() - i - 1 != m.start() + i:
                score[m.end() - i - 1] += coef


def scoring(string, file):
    muscle_cline = MuscleCommandline(input=file)
    stdout, stderr = muscle_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    length = align.get_alignment_length()
    score = [0] * length
    # TODO : manage when the align score is too low
    for record in align:
        pattern = r"%s" % string
        seq = str(record.seq)
        tmp = find_pattern(pattern, seq)
        for m in tmp:
            fill_score_table(score, m, align)
    return score


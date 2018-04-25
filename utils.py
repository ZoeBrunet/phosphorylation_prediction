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
import os
from io import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


def common_parse(parser):
    parser.add_argument('pattern',
                        help='Input Python regular expression you want to detect')
    parser.add_argument('file',
                        help='Input file containing examples')
    parser.add_argument('max_window', type=int, nargs='?', default=15,
                        help='Size of the windows which contain your pattern')


def parse_args_scoring(args):
    parser = argparse.ArgumentParser(description='Run scoring to detect '
                                                 'the frequency of some pattern')
    common_parse(parser)
    return parser.parse_args(args)


def find_pattern(pattern, seq):
    p = re.compile(pattern)
    result = []
    for m in p.finditer(seq):
        result.append(m)
    return result


def fill_score_table(score, m, align, max_window):
    window_length = m.end() - m.start()
    upper = round((window_length / 2))
    if window_length > 0 & window_length <= max_window:
        for i in range(0, upper + 1):
            coef = (i + 1) / (align.__len__() * window_length)
            if m.start() + i <= m.end() - i - 1:
                score[m.start() + i] += coef
            if m.end() - i - 1 > m.start() + i:
                score[m.end() - i - 1] += coef


def run_muscle(file, cluster_name, path):
    muscle_cline = MuscleCommandline(input=file)
    stdout, stderr = muscle_cline()
    file_name = "%s_align.fasta" % cluster_name
    create_align_file = "if [ ! -d %s/align ] ; " \
                        "then mkdir %s/align; " \
                        "fi" % (path, path)
    os.system(create_align_file)
    path2file = "%s/align/%s" %(path, file_name)
    if not os.path.exists(path2file):
        with open(path2file, "w") as g:
            g.write(str(AlignIO.read(StringIO(stdout), "fasta")))
    return AlignIO.read(StringIO(stdout), "fasta")


def relative_position(seq, position):
    j = 0
    for i, char in enumerate(seq):
        if char.isalnum():
            j += 1
        if j >= position:
            break
    return i


def create_window(pos, max_window, length):
    if max_window % 2 == 0:
        max_window += 1
    half_window = (max_window - 1) / 2
    if max_window >= length:
        return [0, length]
    if pos - half_window < 0:
        return [0, max_window]
    if pos + half_window > length:
        return [int(length - max_window), int(length)]
    return [int(pos - half_window), int(pos + half_window)]


def scoring(string, file, max_window):
    align = run_muscle(file)
    length = align.get_alignment_length()
    score = [0] * length
    for record in align:
        pattern = r"%s" % string
        print(pattern)
        seq = str(record.seq)
        tmp = find_pattern(pattern, seq)
        for m in tmp:
            fill_score_table(score, m, align, max_window)
    return score


def score_in_window(file, gene, max_window,
                    pattern, cluster_name, path):
    align = run_muscle(file, cluster_name, path)
    length = align.get_alignment_length()
    pos = None
    for record in align:
        if len(find_pattern(str(gene._get_taxID()), str(record.id))):
            pos = relative_position(record.seq, gene._get_position())
            break
    if pos is not None:
        window = create_window(pos, max_window, length)
        score = [0] * max_window
        for record in align:
            seq = str(record.seq[window[0]:window[1] + 1])
            tmp = find_pattern(pattern, seq)
            for m in tmp:
                fill_score_table(score, m, align, max_window)
    return score


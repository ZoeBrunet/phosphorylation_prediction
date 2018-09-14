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
import math
from source.utils.tools import *


def get_freq_of_pattern(pattern, window, file, max_window):
    score = ["nan"] * max_window
    if os.path.exists(file):
        if len(window):
            with open(file) as f:
                max_window = int(window[2][1] - window[2][0] + 1)
                align = AlignIO.read(f, "fasta")
                score = [0] * max_window
                for record in align:
                    seq = str(record.seq[window[2][0]:window[2][1] + 1])
                    tmp = find_pattern(pattern, seq)
                    for m in tmp:
                        fill_score_table(score, m, align, max_window)
    return score


def get_information_content(window, file):
    info_content = ["nan"] * 3
    if os.path.exists(file):
        if len(window):
            for i, w in enumerate(window):
                summary_align = get_align_info(file)
                info_content[i] = (summary_align.information_content(w[0], w[1],
                                                                     log_base=10,
                                                                     chars_to_ignore=['-', 'X']))
    return info_content


def get_shanon_entropy(window, pssm, max_window):
    shanon_list = ["nan"] * max_window
    if pssm is not None:
        if len(window):
            sub_pssm = []
            if pssm.pssm.__len__() - 1 >= window[2][1]:
                sub_pssm = [pssm[index] for index in range(window[2][0], window[2][1] + 1)]
            else:
                if pssm.pssm.__len__() - 1 >= window[2][0]:
                    sub_pssm = [pssm[index] for index in range(window[2][0], pssm.pssm.__len__())]
            shanon = lambda f: -(f * math.log(f, 2)) if f != 0 else 0
            for i, row in enumerate(sub_pssm):
                tot = sum(row.values())
                shanon_value = 0
                if tot != 0:
                    for val in row.values():
                        freq = val / tot
                        shanon_value += shanon(freq)
                    shanon_list[i] = shanon_value
    return shanon_list


def get_ACH(window, sequence):
    ACH_list = ["nan"] * 3
    if len(window):
        for i, w in enumerate(window):
            ACH = 0.
            seq = str(sequence[int(w[0]):int(w[1] + 1)])
            if "X" in seq:
                ACH_list[i] = "nan"
            else:
                for char in seq:
                    ACH = round(ACH + hydrophobicity[char], 2)
                ACH_list[i] = ACH
    return ACH_list


def get_alignment_ACH(window, file):
    score = ["nan"] * 3
    if os.path.exists(file):
        if len(window):
            score = [0] * 3
            with open(file) as f:
                align = AlignIO.read(f, "fasta")
            for record in align:
                ACH = get_ACH(window, record.seq)
                for i, a in enumerate(ACH):
                    if a != "nan":
                        score[i] += a
            score = [s/align.__len__() for s in score]
    return score

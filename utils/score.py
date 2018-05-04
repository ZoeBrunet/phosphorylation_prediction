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
from utils.tools import *
import math


def get_freq_of_pattern(pattern, window, file):
    with open(file) as f:
        max_window = int(window[1] - window[0] + 1)
        align = AlignIO.read(f, "fasta")
        if len(window):
            score = [0] * max_window
            for record in align:
                seq = str(record.seq[window[0]:window[1] + 1])
                tmp = find_pattern(pattern, seq)
                for m in tmp:
                    fill_score_table(score, m, align, max_window)
    return {"score": score, "nb_align": align.__len__()}


def get_information_content(window, file):
    summary_align = get_align_info(file)
    info_content = summary_align.information_content(window[0], window[1],
                                                     log_base=10,
                                                     chars_to_ignore=['-','X'])
    return info_content


def get_shanon_entropy(window, pssm):
    shanon_list = []
    sub_pssm = [pssm[index] for index in range(window[0], window[1])]
    for row in sub_pssm:
        tot = sum(row.values())
        shanon_value = 0
        shanon = lambda f: -(f * math.log(f, 2)) if f != 0 else 0
        for val in row.values():
            freq = val / tot
            shanon_value += shanon(freq)
        shanon_list.append(shanon_value)
    return shanon_list


def get_ACH(window, sequence):
    ACH = 0.
    if len(window):
        seq = str(sequence[window[0]:window[1] + 1])
        for char in seq:
            ACH = round(ACH + hydrophobicity[char], 2)
    return ACH
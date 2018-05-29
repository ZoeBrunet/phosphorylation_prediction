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
from difflib import SequenceMatcher


def relative_position(seq, position):
    j = -1
    for i, char in enumerate(seq):
        if char.isalnum():
            j += 1
        if j >= position:
            break
    return i


def find_pos_in_alignment(align, sequence, taxID, position):
    pos = None
    for record in align:
        if len(find_pattern(str(taxID), str(record.id))):
            seq = record.seq
            if sequence == str(seq).replace('-', ''):
                pos = relative_position(record.seq, position)
            else:
                match = SequenceMatcher(lambda x: x == "-", sequence,
                                        record.seq).find_longest_match(0,
                                                                       len(sequence),
                                                                       0,
                                                                       len(record.seq))
                if match.a <= position <= match.a + match.size:
                    new_pos = position - match.a
                    pos = match.b + relative_position(record.seq[match.b: match.b + match.size], new_pos)
            break
    return {"position": pos, "sequence": seq}


def create_window(sequence, pos, max_window):
    if pos is not None:
        if max_window % 2 == 0:
            max_window += 1
        half_window = (max_window - 1) / 2
        length = len(sequence)
        if max_window >= length:
            return [[0, pos - 1], [pos + 1, length], [0, length]]
        if pos - half_window < 0:
            return [[0, pos - 1], [pos + 1, max_window - 1], [0, max_window - 1]]
        if pos + half_window > length:
            return [[int(length - max_window) + 1, pos - 1], [pos +1, int(length)],
                    [int(length - max_window) + 1, int(length)]]
        return [[int(pos - half_window), pos - 1], [pos + 1, int(pos + half_window)],
                [int(pos - half_window), int(pos + half_window)]]
    return []

def get_big_window(file):
    with open(file) as f:
        align = AlignIO.read(f, "fasta")
        length = align.get_alignment_length()
    return [0, length]
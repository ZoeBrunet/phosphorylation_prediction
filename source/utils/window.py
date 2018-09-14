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
from source.utils.tools import *
from difflib import SequenceMatcher


def relative_position(seq, position):
    j = -1
    for i, char in enumerate(seq):
        if char.isalnum():
            j += 1
        if j >= position:
            break
    return i


def find_pos_in_alignment(align, sequence, taxID, position, phospho_ELM):
    pos = None
    seq = None
    for record in align:
        if len(find_pattern(str(taxID), str(record.id))):
            seq = record.seq
            rel_seq = str(record.seq).replace('-', '')
            if sequence == rel_seq and phospho_ELM:
                pos = relative_position(record.seq, position)
            else:
                (seq1, seq2) = (rel_seq, sequence) if (len(sequence) < len(rel_seq)) \
                    else (sequence, rel_seq)
                match = SequenceMatcher(None, seq1,
                                        seq2).find_longest_match(0,
                                                                 len(seq1),
                                                                 0,
                                                                 len(seq2))
                (match_align, match_sequence) = (match.a, match.b) if \
                    (len(sequence) < len(rel_seq)) else (match.b, match.a)
                if match_sequence <= position <= match_sequence + match.size - 1 and phospho_ELM:
                    start = relative_position(seq, match_align)
                    new_pos = position - match_sequence
                    pos = start + relative_position(seq[start:], new_pos)
                if not phospho_ELM:
                    if match.size == 13:
                        pos = relative_position(seq, 6 + match_align)
                    else:
                        pos = 6
                        seq = sequence
            break
    return {"position": pos, "sequence": seq}


def create_window(sequence, pos, max_window, phospho_ELM):
    if pos is not None:
        if max_window % 2 == 0:
            max_window += 1
        half_window = (max_window - 1) / 2
        if phospho_ELM:
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
    return [[0, 1], [0, 1], [0, length]]

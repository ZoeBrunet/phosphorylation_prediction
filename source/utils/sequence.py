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
from difflib import SequenceMatcher
from Bio import SeqIO
import os
import math


def find_sequence(path2cluster, seq_list, taxID):
    best_seq = None
    former_id_list = []
    if os.path.exists(path2cluster):
        best_tmp = 0
        for record in SeqIO.parse(open(path2cluster), "fasta"):
            taxonomy = str(record.id).split(":")[0]
            if record.id not in former_id_list:
                former_id_list.append(record.id)
                if float(taxonomy) == float(taxID):
                    tmp = 0
                    for seq in seq_list:
                        match = SequenceMatcher(None, record.seq,
                                                seq).find_longest_match(0, len(record.seq),
                                                                        0, len(seq))
                        if match.size == 13:
                            tmp += 1
                    if tmp > 0:
                        if tmp > best_tmp:
                            best_seq = record.seq
                        if tmp == best_tmp:
                            if len(record.seq) > len(best_seq):
                                best_seq = record.seq
    return str(best_seq)


def is_metazoan(taxID, mt):
    if taxID is not None:
        if not math.isnan(taxID) and taxID is not None:
            info = mt.gettaxon(int(taxID))
            if info is not None:
                if "lineage" in info:
                    if 33208 in info["lineage"]:
                        return True
                    else:
                        return False
    return "nan"

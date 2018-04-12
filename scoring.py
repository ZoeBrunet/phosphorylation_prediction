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
# This program is distributed in the hope that it will be useful,
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import sys
from io import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


def scoring(string, file):
    muscle_cline = MuscleCommandline(input=file)
    stdout, stderr = muscle_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    length = align.get_alignment_length()
    score = [0] * length

    # TODO : manage when the align score is too low

    # TODO : manage when string is longer than 1 char

    for record in align:
        for i, char in enumerate(record):
            if char == string:
                score[i] += 1 / align.__len__()
    return score


if len(sys.argv) > 2:
    print(scoring(sys.argv[1], sys.argv[2]))
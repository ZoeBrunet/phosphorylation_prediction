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
from Bio.Align.Applications import MuscleCommandline


def run_muscle(file_input):
    file = os.path.basename(file_input)
    file_output = "%s_align.fasta" % file[:-6]
    path = os.path.dirname(os.path.dirname(file_input))
    path2align = "%s/align" % path
    if not os.path.exists(path2align):
        os.mkdir(path2align)
    path2outfile = "%s/%s" % (path2align, file_output)
    if not os.path.exists(path2outfile):
        muscle_cline = MuscleCommandline(input=file_input, out=path2outfile)
        os.system(str(muscle_cline))
    return(path2outfile)

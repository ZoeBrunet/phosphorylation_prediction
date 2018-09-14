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
import sys
from source.utils.parser import score_parser
from source.utils.score import get_shanon_entropy
from source.utils.align_ortholog import run_muscle
from source.utils.tools import get_pssm, get_align_info
from source.utils.window import get_big_window

args = score_parser(sys.argv[1:], 'Run shanon_entropy to get the shanon entropy '
                                    'of each position in the alignment')
outputfile = run_muscle(args.file)
summary_align = get_align_info(outputfile)
pssm = get_pssm(summary_align)
window = get_big_window(outputfile)

print(get_shanon_entropy(window, pssm, window[2][1]))

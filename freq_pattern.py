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
from utils.parser import common_parser
from utils.score import freq_of_pattern
from utils.align_ortholog import run_muscle

args = common_parser(sys.argv[1:], 'Run freq_pattern to get the '
                                   'frequency of the pattern')
outputfile = run_muscle(args.file)
window = [0, args.max_window]
print(freq_of_pattern(args.pattern, window, outputfile))

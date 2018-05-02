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

import csv
import sys
from utils.parser import common_parser
from utils.import_csv import import_ortholog
from utils.score import *
from utils.align_ortholog import *
from utils.window import create_window


def align_file(string, file, max_window):
    file_name = os.path.basename(file)
    path = os.path.dirname(file)
    pattern = r"%s" % string
    path2fastas = '%s/fastas' % path
    gene_list = import_ortholog(file, pattern)
    length_list_gene = len(gene_list)
    with open('%s/%s_%s_train_table.csv' % (path, file_name[:-4], string),
              'w', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        writer.writerow(('uniprotID', 'geneID', 'code', 'position',
                         'taxID', 'clusterID', 'sequence', 'freq',
                         'IC'))
        for i, gene in enumerate(gene_list):
            clusterID = gene._get_cluster()
            freq = None
            IC = None
            if clusterID is not None:
                input = "%s.fasta" % gene._get_cluster()
                path2input = '%s/%s' % (path2fastas, input)
                path2align = run_muscle(path2input)
                print_trace(i, length_list_gene, "run muscle")
                with open(path2align) as f:
                    align = AlignIO.read(f, "fasta")
                    length_list_alignment = align.get_alignment_length()
                    window = create_window(max_window, length_list_alignment, align, gene)
                    if len(window):
                        freq = freq_of_pattern(pattern, window, path2align)
                        IC = get_information_content(window, path2align)
            writer.writerow((gene._get_uniprotID(),gene._get_geneID(),
                             gene._get_code(), gene._get_position(),
                             gene._get_taxID(), clusterID,
                             gene._get_sequence(), freq, IC))


args = common_parser(sys.argv[1:], 'Run align to get'
                                   'csv training set')
align_file(args.pattern, args.file, args.max_window)

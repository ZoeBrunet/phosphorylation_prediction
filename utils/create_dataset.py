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
from utils.import_csv import import_ortholog
from utils.score import *
from utils.align_ortholog import *
from utils.window import create_window


def create_training_set(string, file, max_window, phospho_site):
    file_name = os.path.basename(file)
    path = "%s/data" % os.path.abspath(os.path.dirname
                                       (os.path.dirname(__file__)))
    pattern = r"%s" % string
    path2fastas = '%s/fastas' % path
    path2csv = '%s/csv/%s' % (path, string)
    gene_list = import_ortholog(file, pattern, phospho_site)
    length_list_gene = len(gene_list)
    if not os.path.exists(os.path.dirname(path2csv)):
        os.mkdir(os.path.dirname(path2csv))
    if not os.path.exists(path2csv):
        os.mkdir(path2csv)
    suffix = "pos_sites" if phospho_site else "neg_sites"
    with open('%s/%s_%s_%s.csv' % (path2csv, file_name[:-4], string, suffix),
              'w', newline='') as g:
        print(path2csv)
        writer = csv.writer(g, delimiter=";")
        writer.writerow(('uniprotID', 'geneID', 'code', 'position',
                         'taxID', 'clusterID', 'sequence', 'freq',
                         'IC', 'nb_orthologs', 'phosphorylation_site'))
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
                    window = create_window(max_window, length_list_alignment,
                                           align, gene)
                    if len(window):
                        freq = freq_of_pattern(pattern, window, path2align)
                        freq_score = freq["score"]
                        nb_orthologs = freq['nb_align']
                        IC = get_information_content(window, path2align)
            writer.writerow((gene._get_uniprotID(),gene._get_geneID(),
                             gene._get_code(), gene._get_position(),
                             gene._get_taxID(), clusterID,
                             gene._get_sequence(), freq_score, IC,
                             nb_orthologs, gene._get_phosphorylation_site()))

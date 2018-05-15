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


def align_fastas(gene_list, path2fastas, path2align):
    pssm_list = {}
    for i, gene in enumerate(gene_list):
        clusterID = gene._get_cluster()
        length_list_gene = len(gene_list)

        # get align cluster

        if clusterID is not None:
            input = "%s.fasta" % gene._get_cluster()
            path2input = '%s/%s' % (path2fastas, input)
            path2cluster = "%s/%s_align.fasta" % (path2align, input[:-6])

            # align fasta

            if not os.path.exists(path2cluster):
                run_muscle(path2input)
            if clusterID not in pssm_list:
                summary_align = get_align_info(path2cluster)
                pssm_list[clusterID] = get_pssm(summary_align)
        print_trace(i, length_list_gene, "run muscle")
    return pssm_list


def fill_file(gene_list, path2csv, file_name, pattern, string,
              suffix, path2align, max_window, pssm_list):

    # Creation of the output file

    if not os.path.exists(os.path.dirname(path2csv)):
        os.mkdir(os.path.dirname(path2csv))
    if not os.path.exists(path2csv):
        os.mkdir(path2csv)

    # Writing into the output file

    with open('%s/%s_%s_%s.csv' % (path2csv, file_name[:-4], string, suffix),
              'w', newline='') as g:

        # Header

        writer = csv.writer(g, delimiter=";")
        writer.writerow(('uniprotID', 'geneID', 'code', 'position',
                         'taxID', 'clusterID', 'sequence', 'freq_left', 'freq_right', 'freq_tot',
                         'IC_left', 'IC_right', 'IC_tot', 'nb_orthologs', 'phosphorylation_site',
                         'shanon_entropy_left', 'shanon_entropy_right', 'shanon_entropy_tot',
                         'ACH_left', 'ACH_right', 'ACH_tot'))
        length = len(gene_list)
        for i, gene in enumerate(gene_list):
            clusterID = gene._get_cluster()
            freq_score = []
            IC = []
            ACH = []
            shanon_entropy = []
            print_trace(i, length, "fill %s_%s_%s.csv " % (file_name[:-4], string, suffix))
        # get align cluster
            if clusterID is not None:
                input = "%s.fasta" % gene._get_cluster()
                path2cluster = "%s/%s_align.fasta" % (path2align, input[:-6])

                with open(path2cluster) as f:
                    alpha = Alphabet.Gapped(IUPAC.protein)
                    align = AlignIO.read(f, "fasta", alphabet=alpha)
                    length_list_alignment = align.get_alignment_length()
                    window = create_window(max_window, length_list_alignment,
                                           align, gene)
                    nb_orthologs = align.__len__()

                    # Get scores
                if len(window):
                    for w in window:
                            freq_score.append(get_freq_of_pattern(pattern, w, path2cluster))
                            IC.append(get_information_content(w, path2cluster))
                            shanon_entropy.append(get_shanon_entropy(w, pssm_list[clusterID]))
                            ACH.append(get_ACH(w, gene._get_sequence()))

                        # write row

                    writer.writerow((gene._get_uniprotID(), gene._get_geneID(),
                                     gene._get_code(), gene._get_position(),
                                     gene._get_taxID(), clusterID,
                                     gene._get_sequence(), freq_score[0], freq_score[1], freq_score[2],
                                     IC[0], IC[1], IC[2], nb_orthologs, gene._get_phosphorylation_site(),
                                     shanon_entropy[0], shanon_entropy[1], shanon_entropy[2],
                                     ACH[0], ACH[1], ACH[2]))


def create_training_set(string, file, max_window):

    # Initialisation

    file_name = os.path.basename(file)
    path = "%s/data" % os.path.abspath(os.path.dirname
                                       (os.path.dirname(__file__)))
    pattern = r"%s" % string
    path2fastas = '%s/fastas' % path
    path2align = '%s/align' % path
    path2csv = '%s/csv/%s' % (path, string)
    suffix_pos = "pos_sites"
    suffix_neg = "neg_sites"

    # Data importation

    g = import_ortholog(file, pattern)
    gene_list_pos = g["positif"]
    gene_list_neg = g["negatif"]
    pssm_list = align_fastas(gene_list_pos, path2fastas, path2align)
    for gene_list, suffix in zip([gene_list_pos, gene_list_neg], [suffix_pos, suffix_neg]):
        fill_file(gene_list, path2csv, file_name, pattern, string,
                  suffix, path2align, max_window, pssm_list)

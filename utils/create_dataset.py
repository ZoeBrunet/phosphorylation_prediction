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


def print_info(gene, clusterID, freq_score, IC, nb_orthologs, shanon_entropy, ACH,
               writer, window, align):
    writer.writerow(([gene._get_uniprotID(), gene._get_geneID(),
                      gene._get_code(), gene._get_position(),
                      gene._get_taxID(), clusterID,
                      gene._get_sequence(), nb_orthologs, gene._get_phosphorylation_site(),
                      ACH[0], ACH[1], ACH[2], IC[0], IC[1], IC[2], is_metazoan(gene._get_taxID())]
                     + freq_score[2] + shanon_entropy[2]))


    for record in align:
        if len(find_pattern(str(gene._get_taxID()), str(record.id))):
            seq = record.seq
            break
    seq_left = '    '.join(seq[window[0][0]: window[0][1] + 1]) + "    "
    phospho_site = seq[window[0][1] + 1: window[1][0]]
    seq_right = '    '.join(seq[window[1][0]: window[1][1] + 1]) + "    "
    print("\n\033[31;4mInfo\033[0m :")
    print("\n\033[;4mUniprotID\033[0m : %s   \033[;4mGeneID\033[0m : %s   "
          "\033[;4mTaxID\033[0m : %s   \033[;4mPosition\033[0m : %s"
          "   \033[;4mMetazoa\033[0m : %s" % (gene._get_uniprotID(), gene._get_geneID(),
                                               gene._get_taxID(), gene._get_position(),
                                               is_metazoan(gene._get_taxID())))
    print("\nsequence            : \033[34m %s\033[0m%s    \033[32m%s\033[0m \n " % (seq_left, phospho_site, seq_right))

    lamb = lambda n: "NAN" if n == "NA" else round(float(n), 1)
    freq_left = [lamb(element)for element in freq_score[0]]
    freq_phospho = lamb(freq_score[2][len(freq_score[0])])
    freq_right = [lamb(element) for element in freq_score[1]]
    print("freq                : \033[34m %s\033[0m, %s, \033[32m%s\033[0m \n " % (str(freq_left)[1:-1],
                                                str(freq_phospho),
                                                str(freq_right)[1:-1]))
    se_left = [lamb(element) for element in shanon_entropy[0]]
    se_phospho = lamb(shanon_entropy[2][len(shanon_entropy[0])])
    se_right = [lamb(element) for element in shanon_entropy[1]]
    print("shanon entropy      : \033[34m %s\033[0m, %s, \033[32m%s\033[0m \n " % (str(se_left)[1:-1],
                                                                        str(se_phospho),
                                                                        str(se_right)[1:-1]))
    print("information content : \033[34m   %s\033[0m,       %s,       \033[32m%s\033[0m \n " % (IC[0],
                                                                              IC[2],
                                                                              IC[1]))
    print("ACH                 : \033[34m              %s\033[0m,                %s,            "
          "  \033[32m%s\033[0m \n " % (ACH[0], ACH[2], ACH[1]))


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
        header_freq = ["freq_%s" %i for i in range(0, max_window)]
        header_se = ["shanon_entropy_%s" %i for i in range(0, max_window)]
        writer.writerow((['uniprotID', 'geneID', 'code', 'position',
                          'taxID', 'clusterID', 'sequence', 'nb_orthologs', 'phosphorylation_site',
                          'ACH_left', 'ACH_right', 'ACH_tot', 'IC_left', 'IC_right', 'IC_tot', "metazoa"]
                         + header_freq + header_se))
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
                        if len(w):
                            freq_score.append(get_freq_of_pattern(pattern, w, path2cluster))
                            IC.append(get_information_content(w, path2cluster))
                            shanon_entropy.append(get_shanon_entropy(w, pssm_list[clusterID]))
                            ACH.append(get_ACH(w, gene._get_sequence()))

                        # write row

                    print_info(gene, clusterID, freq_score, IC, nb_orthologs, shanon_entropy, ACH,
                               writer, window, align)


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

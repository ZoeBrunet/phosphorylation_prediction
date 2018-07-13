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
import ast
from utils.import_csv import *
from utils.score import *
from utils.align_ortholog import *
from utils.tools import split_fasta
from utils.window import create_window, find_pos_in_alignment


def create_training_set(string, file, max_window, phospho_ELM=True, progression=False,
                        color=False, align_ortho_window=True):
    # Initialisation

    file_name = os.path.basename(file)
    if os.path.basename(os.path.dirname(file)) == "csv" \
            and os.path.basename(os.path.dirname(os.path.dirname(file))) == "data"\
            and os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file))))\
            =="phosphorylation_prediction":
        path = "%s/data" % os.path.abspath(os.path.dirname
                                           (os.path.dirname(__file__)))
    else:
        if not os.path.exists("%s/data" % os.path.dirname(file)):
            os.mkdir("%s/data" % os.path.dirname(file))
        if not os.path.exists("%s/data/csv" % os.path.dirname(file)):
            os.mkdir("%s/data/csv" % os.path.dirname(file))
        if not os.path.exists("%s/data/csv/%s" % (os.path.dirname(file), string)):
            os.mkdir("%s/data/csv/%s" % (os.path.dirname(file), string))
        path = "%s/data" % os.path.dirname(file)

    pattern = r"%s" % string
    path2fastas = '%s/fastas' % path
    path2align = '%s/align' % path
    path2csv = '%s/csv/%s' % (path, string)

    # Data importation

    index_file = import_ortholog(file, pattern, phospho_ELM, progression)
    genes = pd.read_csv(index_file, sep=';')

    # Creation of csv

    with open('%s/%s_%s_phospho_sites.csv' % (path2csv, file_name[:-4], string),
              'w', newline='') as g:
        writer = csv.writer(g, delimiter=";")

        # Header

        header_freq_metazoa = ["freq_metazoa_%s" % i for i in range(0, max_window)]
        header_se_metazoa = ["shanon_entropy_metazoa_%s" % i for i in range(0, max_window)]
        header_freq_non_metazoa = ["freq_metazoa_%s" % i for i in range(0, max_window)]
        header_se_non_metazoa = ["shanon_entropy_metazoa_%s" % i for i in range(0, max_window)]
        writer.writerow((['uniprotID', 'geneID', 'position',
                          'taxID', 'clusterID', 'sequence', 'seq_in_window',
                          'nb_orthologs', 'phosphorylation_site',
                          'ACH_prot_left', 'ACH_prot_right', 'ACH_prot_tot',
                          'ACH_metazoa_left', 'ACH_metazoa_right', 'ACH_metazoa_tot',
                          'ACH_non_metazoa_left', 'ACH_non_metazoa_right', 'ACH_non_metazoa_tot',
                          'IC_metazoa_left', 'IC_metazoa_right', 'IC_metazoa_tot',
                          'IC_non_metazoa_left', 'IC_non_metazoa_right', 'IC_non_metazoa_tot',
                          "metazoa"]
                         + header_freq_metazoa + header_freq_non_metazoa +
                         header_se_metazoa + header_se_non_metazoa))
        length = len(genes)
        gene_seq = genes["sequence"] if phospho_ELM else genes["seq_in_window"]
        gene_pos = zip(genes["pos_sites"], genes["neg_sites"]) if phospho_ELM else genes["pos_sites"]
        for i, (uniprotID, geneID, taxID,
                metazoan, sequence, sites, clusterID) in enumerate(zip(genes["uniprotID"], genes["geneID"],
                                                       genes["taxID"], genes["metazoan"],
                                                       gene_seq, gene_pos, genes["clusterID"])):
            window_seq = None if phospho_ELM else sequence
            input = "%s.fasta" % clusterID
            path2input = '%s/%s' % (path2fastas, input)
            path2cluster = "%s/%s.fasta" % (path2fastas, input[:-6])
            path2aligncluster = "%s/%s_align.fasta" % (path2align, input[:-6])
            path2alignclustermetazoa = "%s/%s_align_metazoa.fasta" % (path2align, input[:-6])
            path2alignclusternonmetazoa = "%s/%s_align_nonmetazoa.fasta" % (path2align, input[:-6])

            # align fasta
            pssm_metazoa = None
            pssm_non_metazoa = None
            ortholog_metazoa = False
            ortholog_non_metazoa = False
            if os.path.exists(path2cluster):
                if not os.path.exists(path2aligncluster):
                    run_muscle(path2input)
                if ((not os.path.exists(path2alignclustermetazoa) or
                    not os.path.exists(path2alignclusternonmetazoa)) and
                        os.path.exists(path2aligncluster)):
                    split_fasta(path2aligncluster)
                if os.path.exists(path2alignclustermetazoa):
                    summary_align_metazoa = get_align_info(path2alignclustermetazoa)
                    pssm_metazoa = get_pssm(summary_align_metazoa)
                    ortholog_metazoa = True
                if os.path.exists(path2alignclusternonmetazoa):
                    summary_align_non_metazoa = get_align_info(path2alignclusternonmetazoa)
                    pssm_non_metazoa = get_pssm(summary_align_non_metazoa)
                    ortholog_non_metazoa = True

            # fill csv
            phospho_bool = [True, False] if phospho_ELM else [True]
            if not phospho_ELM:
                sites = [sites]
            for position_list, phosphorylation_site in zip(sites, phospho_bool):
                for position in ast.literal_eval(position_list):
                    if progression:
                        print_trace(i, length, "filling CSV file ")
                    window = create_window(sequence, position, max_window, phospho_ELM)
                    rel_window = []
                    nb_orthologs = 0
                    if not phospho_ELM:
                        sequence = None
                    rel_sequence = sequence
                    #for vecteur bool pour metazoa et non metazoa zip
                    if ortholog_metazoa or ortholog_non_metazoa:
                        with open(path2aligncluster) as f:
                            alpha = Alphabet.Gapped(IUPAC.protein)
                            align = AlignIO.read(f, "fasta", alphabet=alpha)
                            nb_orthologs = align.__len__()
                            if not phospho_ELM:
                                sequence = find_seq(align, taxID)
                            finder = find_pos_in_alignment(align, sequence, taxID,
                                                           position, phospho_ELM) if phospho_ELM \
                                else find_pos_in_alignment(align, window_seq, taxID,
                                                           position, phospho_ELM)
                            rel_pos = finder["position"]
                            rel_sequence = finder["sequence"]
                            rel_window = create_window(rel_sequence, rel_pos, max_window, True)

                    # Score orthologs

                    freq_metazoa = get_freq_of_pattern(pattern, rel_window, path2alignclustermetazoa, max_window)
                    freq_non_metazoa = get_freq_of_pattern(pattern, rel_window, path2alignclusternonmetazoa, max_window)
                    shanon_entropy_metazoa = get_shanon_entropy(rel_window, pssm_metazoa, max_window)
                    shanon_entropy_non_metazoa = get_shanon_entropy(rel_window, pssm_non_metazoa, max_window)
                    IC_metazoa = get_information_content(rel_window, path2alignclustermetazoa)
                    IC_non_metazoa = get_information_content(rel_window, path2alignclusternonmetazoa)

                    # Score sequence

                    ACH_prot = get_ACH(window, sequence) if sequence is not None \
                        else get_ACH(create_window(window_seq, 6, 13, phospho_ELM), window_seq)
                    ACH_metazoa = get_alignment_ACH(rel_window, path2alignclustermetazoa)
                    ACH_non_metazoa = get_alignment_ACH(rel_window, path2alignclusternonmetazoa)


                    # Fill csv
                    if window_seq is None:
                        window_seq = sequence[window[0][0]: window[1][1] + 1]
                    writer.writerow([uniprotID, geneID, position, taxID, clusterID,
                                     sequence, window_seq, nb_orthologs,
                                     phosphorylation_site, ACH_prot[0], ACH_prot[1],
                                     ACH_prot[2], ACH_metazoa[0], ACH_metazoa[1],
                                     ACH_metazoa[2], ACH_non_metazoa[0], ACH_non_metazoa[1],
                                     ACH_non_metazoa[2], IC_metazoa[0], IC_metazoa[1], IC_metazoa[2],
                                     IC_non_metazoa[0], IC_non_metazoa[1], IC_non_metazoa[2],
                                     metazoan] + freq_metazoa + freq_non_metazoa +
                                    shanon_entropy_metazoa + shanon_entropy_non_metazoa)

                    # Print infos

                    if not len(rel_window):
                        rel_window = window
                        rel_sequence = sequence
                    if sequence is None:
                        rel_window = create_window(window_seq, 6, 13, phospho_ELM)
                        rel_sequence = window_seq

                    red = "\n\033[31;4m" if color else "\n"
                    blue = "\033[34m" if color else ""
                    green = "\033[32m" if color else ""
                    white = "\033[37m" if color else ""
                    underline = "\033[;4m" if color else ""
                    end = "\033[0m" if color else ""
                    print("%sInfo%s :" % (red, end))
                    print("\n%sUniprotID%s : %s   %sGeneID%s : %s   "
                          "%sTaxID%s : %s   %sPosition%s : %s"
                          "   %sMetazoa%s : %s   %sPhosphorylation%s : %s"
                          % (underline, end, uniprotID, underline, end, geneID,
                             underline, end, taxID, underline, end, position,
                             underline, end, metazoan, underline, end, phosphorylation_site))

                    half_window = int((max_window - 1) / 2)
                    space = [20, 5*half_window, 5, 5*half_window]

                    # Print sequence
                    for name, freq, shanon_entropy, IC, ACH in zip(["metazoa", "non metazoa"],
                                                                   [freq_metazoa, freq_non_metazoa],
                                                                   [shanon_entropy_metazoa, shanon_entropy_non_metazoa],
                                                                   [IC_metazoa, IC_non_metazoa],
                                                                   [ACH_metazoa, ACH_non_metazoa]):
                        print("%s%s%s :" % (red, name, end))
                        seq_left = (' ' * 4).join(rel_sequence[rel_window[0][0]:
                                                               rel_window[0][1] + 1])+(" " * 4)
                        phospho_site = (' ' * 4).join(rel_sequence[rel_window[0][1] + 1:
                                                                   rel_window[1][0]])+(" " * 4)
                        seq_right = (' ' * 4).join(rel_sequence[rel_window[1][0]:
                                                                rel_window[1][1] + 1]) + (" " * 4)
                        print("\nsequence%s:%s%s%s%s%s%s%s \n "
                              % (" " * (space[0] - len("sequence")), blue, seq_left, end,
                                 phospho_site, green, seq_right, end))

                        if align_ortho_window:
                            align = AlignIO.read(path2aligncluster, "fasta")
                            for record in align:
                                if not len(find_pattern(str(taxID), str(record.id))):
                                    seq_left = (' ' * 4).join(record.seq[rel_window[0][0]:
                                                                         rel_window[0][1] + 1]) + (" " * 4)
                                    phospho_site = (' ' * 4).join(record.seq[rel_window[0][1] + 1:
                                                                               rel_window[1][0]]) + (" " * 4)
                                    seq_right = (' ' * 4).join(record.seq[rel_window[1][0]:
                                                                            rel_window[1][1] + 1]) + (" " * 4)
                                    print("%s%s:%s%s%s%s%s%s%s \n "
                                          % (record.id, " " * (space[0] - len(record.id)), blue, seq_left, end,
                                             phospho_site, green, seq_right, end))

                        # Print frequence

                        lamb = lambda n: "nan" if n == "nan" else round(float(n), 1)
                        freq_left = [lamb(element) for element in freq[0: rel_window[0][1] - rel_window[0][0] + 1]]
                        freq_phospho = lamb(freq[rel_window[0][1] - rel_window[0][0] + 1])
                        freq_right = [lamb(element) for element in freq[rel_window[0][1] -
                                                                            rel_window[0][0] + 2: max_window]]
                        print("%sfreq%s:%s%s%s%s, %s ,%s%s%s "
                              % (white, " " * (space[0] - len("freq")), end, blue,
                                 str(freq_left)[1:-1], end, str(freq_phospho), green,
                                 str(freq_right)[1:-1], end))

                        # Print shanon entropy

                        se_left = [lamb(element) for element in shanon_entropy[0:
                                                                               rel_window[0][1]
                                                                               - rel_window[0][0] + 1]]
                        se_phospho = lamb(shanon_entropy[rel_window[0][1] - rel_window[0][0] + 1])
                        se_right = [lamb(element) for element in shanon_entropy[rel_window[0][1] -
                                                                                rel_window[0][0] + 2:
                                                                                max_window]]
                        print("%sshanon entropy%s:%s%s%s%s, %s, %s%s%s \n "
                              % (white, " " * (space[0] - len("shanon entropy")), end, blue, str(se_left)[1:-1],
                                 end, str(se_phospho), green, str(se_right)[1:-1], end))

                        print("%sinformation content :%s%s%s%s%s%s,%s,%s%s%s%s\n "
                              % (white, " " * (int(space[1] / 2) - 3), end, blue,
                                 lamb(IC[0]), " " * (int(space[1] / 2) - 3), end,
                                 lamb(IC[2]), " " * (int(space[1] / 2) - 3),
                                 green, lamb(IC[1]), end))
                        print("ACH alignment%s:%s%s%s%s%s,%s,%s%s%s%s \n "
                              % (" " * (space[0] - len("ACH alignment")), blue,
                                 " " * (int(space[1] / 2) - 3), lamb(ACH[0]), end,
                                 " " * (int(space[1] / 2) - 3),
                                 lamb(ACH[2]), green,
                                 " " * (int(space[1] / 2) - 3), lamb(ACH[1]), end))
                        print("ACH prot%s:%s%s%s%s%s,%s,%s%s%s%s \n "
                              % (" " * (space[0] - len("ACH prot")), blue,
                                 " " * (int(space[1] / 2) - 3), lamb(ACH_prot[0]), end,
                                 " " * (int(space[1] / 2) - 3),
                                 lamb(ACH_prot[2]), green,
                                 " " * (int(space[1] / 2) - 3), lamb(ACH_prot[1]), end))
        if not phospho_ELM:
            df = pd.read_csv("%s/table_%s_phospho_sites.csv" % (path2csv, string), sep=';')
            for index, row in df.iterrows():
                if not row["phosphorylation_site"]:
                    writer.writerow(row)

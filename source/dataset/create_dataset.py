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
import queue
import threading
import pandas as pd
import csv
import numpy as np
from threading import RLock
from source.utils.score import *
from source.utils.window import create_window, find_pos_in_alignment


class Gene:

    def __init__(self, uniprotID, geneID, position, taxID,
                 clusterID, sequence, window_seq, nb_orthologs,
                 metazoan, nb_orthologs_metazoa, nb_orthologs_non_metazoa,
                 pssm_non_metazoa, pssm_metazoa, pattern):

        self.uniprotID = uniprotID
        self.geneID = geneID
        self.positions = [position]
        self.taxID = taxID
        self.clusterID = clusterID
        self.sequence = sequence
        self.window_seq = window_seq
        self.nb_orthologs = nb_orthologs
        self.metazoan = metazoan
        self.nb_orthologs_metazoa = nb_orthologs_metazoa
        self.nb_orthologs_non_metazoa = nb_orthologs_non_metazoa
        self.pssm_non_metazoa = pssm_non_metazoa
        self.pssm_metazoa = pssm_metazoa
        self.pattern = set(pattern)

    def update_positions(self, new_position):
        self.positions.append(new_position)

    def update_pattern(self, new_pattern):
        self.pattern.add(new_pattern)

    def _get_uniprotID(self):
        return self.uniprotID

    def _get_geneID(self):
        return self.geneID

    def _get_positions(self):
        return self.positions

    def _get_taxID(self):
        return self.taxID

    def _get_clusterID(self):
        return self.clusterID

    def _get_sequence(self):
        return self.sequence

    def _get_window_seq(self):
        return self.window_seq

    def _get_nb_orthologs(self):
        return self.nb_orthologs

    def _get_metazoan(self):
        return self.metazoan

    def _get_nb_orthologs_metazoa(self):
        return self.nb_orthologs_metazoa

    def _get_nb_orthologs_non_metazoa(self):
        return self.nb_orthologs_non_metazoa

    def _get_pssm_non_metazoa(self):
        return self.pssm_non_metazoa

    def _get_pssm_metazoa(self):
        return self.pssm_metazoa

    def _get_pattern(self):
        return self.pattern


lock = RLock()


class fill_table(Thread):
    def __init__(self, genes, phospho_ELM, path2fastas, path2align, max_window,
                 path, pattern, color, align_ortho_window, writer):
        Thread.__init__(self)
        self.genes = genes
        self.phospho_ELM = phospho_ELM
        self.path2fastas = path2fastas
        self.path2align = path2align
        self.max_window = max_window
        self.path = path
        self.pattern = pattern
        self.color = color
        self.align_ortho_window = align_ortho_window
        self.writer = writer

    def run(self):
        genes = self.genes
        phospho_ELM = self.phospho_ELM
        path2fastas = self.path2fastas
        path2align = self.path2align
        max_window = self.max_window
        path = self.path
        pattern = self.pattern
        color = self.color
        alpha = Alphabet.Gapped(IUPAC.protein)

        gene_seq = genes["sequence"] if phospho_ELM else zip(genes["seq_in_window"], genes["neg_seq_in_window"])
        gene_pos = zip(genes["pos_sites"], genes["neg_sites"])

        for i, (uniprotID, geneID, taxID,
                metazoan, sequence_list, sites, clusterID, seq) in enumerate(zip(genes["uniprotID"], genes["geneID"],
                                                                                 genes["taxID"], genes["metazoan"],
                                                                                 gene_seq, gene_pos, genes["clusterID"],
                                                                                 genes["sequence"])):

            nb_orthologs_metazoa = 0
            nb_orthologs_non_metazoa = 0
            path2cluster = "%s/%s.fasta" % (path2fastas, clusterID)
            path2aligncluster = "%s/%s_align.fasta" % (path2align, clusterID)
            path2alignclustermetazoa = "%s/metazoa/%s_align_metazoa.fasta" % (path, clusterID)
            path2alignclusternonmetazoa = "%s/non_metazoa/%s_align_non_metazoa.fasta" % (path, clusterID)

            # align fasta
            pssm_metazoa = None
            pssm_non_metazoa = None
            ortholog_metazoa = False
            ortholog_non_metazoa = False
            if os.path.exists(path2cluster):
                taxid_list = []
                os.makedirs("%s/sorted_fastas" % path, exist_ok=True)
                sorted_file = "%s/sorted_fastas/%s_sorted.fasta" % (path, clusterID)
                seq_find = False
                sorted = os.path.exists(sorted_file)
                aligned = os.path.exists(path2aligncluster)
                split = os.path.exists(path2alignclustermetazoa) or os.path.exists(path2alignclusternonmetazoa)
                if not sorted:
                    for record in SeqIO.parse(open(path2cluster), "fasta"):
                        to_add = True
                        taxonomy = str(record.id).split(":")[0]
                        if float(taxonomy) == float(taxID) and taxID not in taxid_list and not seq_find:
                            match = SequenceMatcher(None, record.seq,
                                                    sequence).find_longest_match(0, len(record.seq),
                                                                                 0, len(sequence))
                            if match.size == 13:
                                add_seq(record, sorted_file, taxid_list, taxonomy)
                                seq_find = True
                            else:
                                to_add = False
                        if taxonomy not in taxid_list and to_add:
                            add_seq(record, sorted_file, taxid_list, taxonomy)
                    if not seq_find:
                        os.remove(sorted_file)
                        print("No corresponding sequence in fasta for %s" % clusterID)
                    else:
                        sorted = True
                if sorted and not aligned:
                    t = threading.Thread(target=run_muscle, args=(sorted_file,))
                    t.start()
                    t.join()
                if aligned and not split:
                    split_fasta(path2aligncluster)
                if os.path.exists(path2alignclustermetazoa):
                    with open(path2alignclustermetazoa) as f_metazoa:
                        align_metazoa = AlignIO.read(f_metazoa, "fasta", alphabet=alpha)
                        nb_orthologs_metazoa = align_metazoa.__len__()
                    summary_align_metazoa = get_align_info(path2alignclustermetazoa)
                    pssm_metazoa = get_pssm(summary_align_metazoa)
                    ortholog_metazoa = True
                if os.path.exists(path2alignclusternonmetazoa):
                    with open(path2alignclusternonmetazoa) as f_non_metazoa:
                        align_non_metazoa = AlignIO.read(f_non_metazoa, "fasta", alphabet=alpha)
                        nb_orthologs_non_metazoa = align_non_metazoa.__len__()
                    summary_align_non_metazoa = get_align_info(path2alignclusternonmetazoa)
                    pssm_non_metazoa = get_pssm(summary_align_non_metazoa)
                    ortholog_non_metazoa = True

            # fill csv
            phospho_bool = [True, False]
            for position_list, phosphorylation_site, seq_list in zip(sites, phospho_bool, sequence_list):
                for position, sequence in zip(ast.literal_eval(position_list), ast.literal_eval(seq_list)):
                    seq_csv = seq
                    window_seq = None if phospho_ELM else sequence
                    que = queue.Queue()
                    window = function_in_thread(que, [sequence,
                                                      position,
                                                      max_window,
                                                      phospho_ELM], create_window)
                    rel_window = []
                    nb_orthologs = 0
                    #for vecteur bool pour metazoa et non metazoa zip
                    if ortholog_metazoa or ortholog_non_metazoa:
                        with open(path2aligncluster) as f:
                            align = AlignIO.read(f, "fasta", alphabet=alpha)
                            nb_orthologs = align.__len__()
                            finder = function_in_thread(que, [align, seq, taxID,
                                                              position, True], find_pos_in_alignment)
                            rel_position = finder["position"]
                            rel_sequence = finder["sequence"]
                            if rel_sequence != '':
                                rel_window = function_in_thread(que, [rel_sequence, rel_position, max_window, True],
                                                                create_window)
                            else:
                                rel_sequence = sequence
                                seq_csv = None
                    # Score orthologs
                    que = queue.Queue()

                    freq_metazoa = function_in_thread(que, [pattern, rel_window, path2alignclustermetazoa, max_window],
                                                      get_freq_of_pattern)
                    freq_non_metazoa = function_in_thread(que, [pattern, rel_window, path2alignclusternonmetazoa, max_window],
                                                          get_freq_of_pattern)
                    shanon_entropy_metazoa = function_in_thread(que, [rel_window, pssm_metazoa, max_window],
                                                                get_shanon_entropy)
                    shanon_entropy_non_metazoa = function_in_thread(que, [rel_window, pssm_non_metazoa, max_window],
                                                                    get_shanon_entropy)
                    IC_metazoa = function_in_thread(que, [rel_window, path2alignclustermetazoa],
                                                    get_information_content)
                    IC_non_metazoa = function_in_thread(que, [rel_window, path2alignclusternonmetazoa],
                                                        get_information_content)

                    # Score sequence

                    ACH_prot = function_in_thread(que, [window, sequence], get_ACH) if sequence is not None \
                        else function_in_thread(que, [create_window(window_seq, 6, 13, phospho_ELM), window_seq],
                                                get_ACH)
                    ACH_metazoa = function_in_thread(que, [rel_window, path2alignclustermetazoa],
                                                     get_alignment_ACH)
                    ACH_non_metazoa = function_in_thread(que, [rel_window, path2alignclusternonmetazoa],
                                                         get_alignment_ACH)

                    with lock:
                    # Fill csv
                        if window_seq is None:
                            window_seq = sequence[window[0][0]: window[1][1] + 1]
                        self.writer.writerow([uniprotID, geneID, position, taxID, clusterID,
                                              seq_csv, window_seq, nb_orthologs, nb_orthologs_metazoa,
                                              nb_orthologs_non_metazoa, phosphorylation_site, ACH_prot[0],
                                              ACH_prot[1], ACH_prot[2], ACH_metazoa[0], ACH_metazoa[1],
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
                            rel_window = function_in_thread(que, [window_seq, 6, 13, phospho_ELM], create_window)
                            rel_sequence = window_seq
                        if rel_sequence=='':
                            rel_sequence = seq_csv

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
                        for name, freq, shanon_entropy, IC, ACH, nb_ortho in zip(["metazoa", "non metazoa"],
                                                                                    [freq_metazoa, freq_non_metazoa],
                                                                                    [shanon_entropy_metazoa,
                                                                                     shanon_entropy_non_metazoa],
                                                                                    [IC_metazoa, IC_non_metazoa],
                                                                                    [ACH_metazoa, ACH_non_metazoa],
                                                                                    [nb_orthologs_metazoa,
                                                                                     nb_orthologs_non_metazoa]):
                            print("%s%s%s :    %snb_ortholog:%s %s" % (red, name, end, underline,
                                                                       end, nb_ortho))
                            seq_left = (' ' * 4).join(rel_sequence[rel_window[0][0]:
                                                                   rel_window[0][1] + 1])+(" " * 4)
                            phospho_site = (' ' * 4).join(rel_sequence[rel_window[0][1] + 1:
                                                                       rel_window[1][0]])+(" " * 4)
                            seq_right = (' ' * 4).join(rel_sequence[rel_window[1][0]:
                                                                    rel_window[1][1] + 1]) + (" " * 4)
                            print("\nsequence%s:%s%s%s%s%s%s%s\n "
                                  % (" " * (space[0] - len("sequence")), blue, seq_left, end,
                                     phospho_site, green, seq_right, end))

                            if self.align_ortho_window:
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
                            print("%sfreq%s:%s%s%s%s, %s, %s%s%s \n"
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
                            print("%sshanon entropy%s:%s%s%s%s, %s, %s%s%s\n "
                                  % (white, " " * (space[0] - len("shanon entropy")), end, blue, str(se_left)[1:-1],
                                     end, str(se_phospho), green, str(se_right)[1:-1], end))

                            print("%sinformation content :%s%s%s%s%s%s,%s,%s%s%s%s\n "
                                  % (white, " " * (int(space[1] / 2) - 3), end, blue,
                                     lamb(IC[0]), " " * (int(space[1] / 2) - 3), end,
                                     lamb(IC[2]), " " * (int(space[1] / 2) - 3),
                                     green, lamb(IC[1]), end))
                            print("ACH alignment%s:%s%s%s%s%s,%s,%s%s%s%s\n "
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


def create_training_set(string, max_window, nthread, file, phospho_ELM=True,
                        color=False, align_ortho_window=True, output_file=None):
    # Initialisation

    if os.path.basename(os.path.dirname(file)) == "%s" % string \
            and os.path.basename(os.path.dirname(os.path.dirname(file))) == "csv" \
            and os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file)))) == "data":
        path = "%s" % os.path.abspath(os.path.dirname(os.path.dirname
                                      (os.path.dirname(file))))
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

    genes = pd.read_csv(file, sep=';')

    # Creation of csv

    csv_dataset = '%s/%s_phospho_sites.csv' % (path2csv, string)
    if output_file is not None:
        csv_dataset = output_file
    if not os.path.exists(csv_dataset):
        os.system('touch %s' % csv_dataset)
    with open(csv_dataset, 'a+', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        g.seek(0)
        first_char = g.read(1)

        # Header
        if not first_char:
            header_freq_metazoa = ["freq_metazoa_%s" % i for i in range(0, max_window)]
            header_se_metazoa = ["shanon_entropy_metazoa_%s" % i for i in range(0, max_window)]
            header_freq_non_metazoa = ["freq_nonmetazoa_%s" % i for i in range(0, max_window)]
            header_se_non_metazoa = ["shanon_entropy_nonmetazoa_%s" % i for i in range(0, max_window)]
            writer.writerow((['uniprotID', 'geneID', 'position',
                              'taxID', 'clusterID', 'sequence', 'seq_in_window',
                              'nb_orthologs', 'nb_orthologs_metazoa', 'nb_orthologs_nonmetazoa',
                              'phosphorylation_site','ACH_prot_left', 'ACH_prot_right',
                              'ACH_prot_tot', 'ACH_metazoa_left', 'ACH_metazoa_right',
                              'ACH_metazoa_tot', 'ACH_nonmetazoa_left', 'ACH_nonmetazoa_right',
                              'ACH_nonmetazoa_tot', 'IC_metazoa_left', 'IC_metazoa_right',
                              'IC_metazoa_tot', 'IC_nonmetazoa_left', 'IC_nonmetazoa_right',
                              'IC_nonmetazoa_tot', 'metazoa']
                             + header_freq_metazoa + header_freq_non_metazoa +
                             header_se_metazoa + header_se_non_metazoa))
        data_thread = np.array_split(genes, nthread)
        thread_list = []
        for data in data_thread:
            thread_list.append(fill_table(data, phospho_ELM,
                                          path2fastas, path2align,
                                          max_window, path, pattern,
                                          color, align_ortho_window,
                                          writer))
        for thread in thread_list:
            thread.start()
        for thread in thread_list:
            thread.join()

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
import threading
from threading import Thread, RLock
import queue
import csv
from utils.window import *
from utils.window import create_window
from utils.score import *

lock = RLock()


def function_in_thread(que, args, func):
    t = Thread(target=lambda q, arg: q.put(func(*arg)), args=(que, args[:]))
    t.start()
    t.join()
    result = None
    while not que.empty():
        result = que.get()
    return result


class fill_neg_sites(Thread):
    def __init__(self, pattern, phospho_sites, path2align, max_window, color, path, writer,
                 align_ortho_window):
        Thread.__init__(self)
        self.pattern = pattern
        self.phospho_sites = phospho_sites
        self.path2align = path2align
        self.max_window = max_window
        self.color = color
        self.path = path
        self.writer = writer
        self.align_ortho_window = align_ortho_window

    def run(self):
        pattern = self.pattern
        phospho_sites = self.phospho_sites
        path2align = self.path2align
        max_window = self.max_window
        color = self.color
        path = self.path
        writer = self.writer
        alpha = Alphabet.Gapped(IUPAC.protein)

        for id, val in phospho_sites.items():
            for position in val._get_positions():
                clusterID = val._get_clusterID()
                window_seq = val._get_window_seq()
                sequence = val._get_sequence()
                taxID = val._get_taxID()
                uniprotID = val._get_uniprotID()
                metazoan = val._get_metazoan()
                geneID = val._get_geneID()
                nb_orthologs = val._get_nb_orthologs()
                nb_orthologs_non_metazoa = val._get_nb_orthologs_non_metazoa()
                nb_orthologs_metazoa = val._get_nb_orthologs_metazoa()
                pssm_metazoa = val._get_pssm_metazoa()
                pssm_non_metazoa = val._get_pssm_non_metazoa()
                p = re.compile(pattern)
                result = []

                for m in p.finditer(sequence):
                    result.append(m)
                for m in result:
                    new_position = round((m.end() + m.start() - 1) / 2)
                    neg = True
                    if abs(new_position - position) <= 50:
                        neg = False
                    if neg:
                        path2aligncluster = "%s/%s_align.fasta" % (path2align, clusterID)
                        path2alignclustermetazoa = "%s/metazoa/%s_align_metazoa.fasta" % (path, clusterID)
                        path2alignclusternonmetazoa = "%s/non_metazoa/%s_align_non_metazoa.fasta" % (path, clusterID)

                        window = create_window(sequence, new_position, max_window, False)
                        rel_window = []
                        ortholog_metazoa = os.path.exists(path2alignclustermetazoa)
                        ortholog_non_metazoa = os.path.exists(path2alignclusternonmetazoa)
                        if ortholog_metazoa or ortholog_non_metazoa:
                            with open(path2aligncluster) as f:
                                align = AlignIO.read(f, "fasta", alphabet=alpha)
                                que = queue.Queue()
                                t = Thread(target=lambda q, align,
                                                         sequence,
                                                         taxID,
                                                         new_position: q.put(find_pos_in_alignment(align,
                                                                                                   sequence,
                                                                                                   taxID,
                                                                                                   new_position,
                                                                                                   True)), args=(que,
                                                                                                                 align,
                                                                                                                 sequence,
                                                                                                                 taxID,
                                                                                                                 new_position))
                                t.start()
                                t.join()
                                finder = None
                                while not que.empty():
                                    finder = que.get()
                                rel_pos = 0 if "position" not in finder else finder["position"]
                                rel_sequence = "" if "position" not in finder else finder["sequence"]
                                rel_window = create_window(rel_sequence, rel_pos, max_window, True)

                                # Score orthologs
                        que = queue.Queue()
                        freq_metazoa = function_in_thread(que, [pattern, rel_window,
                                                                path2alignclustermetazoa, max_window],
                                                          get_freq_of_pattern)
                        freq_non_metazoa = function_in_thread(que, [pattern, rel_window,
                                                                path2alignclusternonmetazoa, max_window],
                                                          get_freq_of_pattern)
                        shanon_entropy_metazoa = function_in_thread(que, [rel_window, pssm_metazoa,
                                                                          max_window],
                                                                    get_shanon_entropy)
                        shanon_entropy_non_metazoa = function_in_thread(que, [rel_window, pssm_non_metazoa,
                                                                              max_window],
                                                                        get_shanon_entropy)
                        IC_metazoa = function_in_thread(que, [rel_window, path2alignclustermetazoa],
                                                        get_information_content)
                        IC_non_metazoa = function_in_thread(que, [rel_window, path2alignclusternonmetazoa],
                                                        get_information_content)

                        # Score sequence

                        ACH_prot = function_in_thread(que, [window, sequence],
                                                      get_ACH)
                        ACH_metazoa = function_in_thread(que, [rel_window, path2alignclustermetazoa],
                                                         get_alignment_ACH)
                        ACH_non_metazoa = function_in_thread(que, [rel_window, path2alignclusternonmetazoa],
                                                         get_alignment_ACH)

                        with lock:
                            # Fill csv
                            writer.writerow([uniprotID, geneID, new_position, taxID, clusterID,
                                             sequence, window_seq, nb_orthologs, nb_orthologs_metazoa,
                                             nb_orthologs_non_metazoa, False, ACH_prot[0],
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

                            red = "\n\033[31;4m" if color else "\n"
                            blue = "\033[34m" if color else ""
                            green = "\033[32m" if color else ""
                            white = "\033[37m" if color else ""
                            underline = "\033[;4m" if color else ""
                            end = "\033[0m" if color else ""
                            print("%sInfo%s :" % (red, end))
                            print("\n%sUniprotID%s : %s   %sGeneID%s : %s   %sTaxID%s : %s   %sPosition%s : %s"
                                  "   %sMetazoa%s : %s   %sPhosphorylation%s : %s"
                                  % (underline, end, uniprotID, underline, end, geneID, underline, end, taxID, underline,
                                     end, new_position, underline, end, metazoan, underline, end, False))

                            half_window = int((max_window - 1) / 2)
                            space = [20, 5 * half_window, 5, 5 * half_window]

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
                                                                       rel_window[0][1] + 1]) + (" " * 4)
                                phospho_site = (' ' * 4).join(rel_sequence[rel_window[0][1] + 1:
                                                                           rel_window[1][0]]) + (" " * 4)
                                seq_right = (' ' * 4).join(rel_sequence[rel_window[1][0]:
                                                                        rel_window[1][1] + 1]) + (" " * 4)
                                print("\nsequence%s:%s%s%s%s%s%s%s\n " % (" " * (space[0] - len("sequence")),
                                                                          blue, seq_left, end, phospho_site,
                                                                          green, seq_right, end))

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
                                freq_left = [lamb(element) for element in
                                             freq[0: rel_window[0][1] - rel_window[0][0] + 1]]
                                freq_phospho = lamb(freq[rel_window[0][1] - rel_window[0][0] + 1])
                                freq_right = [lamb(element) for element in freq[rel_window[0][1] -
                                                                                rel_window[0][0] + 2: max_window]]
                                print("%sfreq%s:%s%s%s%s, %s, %s%s%s \n" % (white, " " * (space[0] - len("freq")),
                                                                            end, blue, str(freq_left)[1:-1], end,
                                                                            str(freq_phospho), green,
                                                                            str(freq_right)[1:-1], end))

                                # Print shanon entropy

                                se_left = [lamb(element) for element in shanon_entropy[0:rel_window[0][1]
                                                                                         - rel_window[0][0] + 1]]
                                se_phospho = lamb(shanon_entropy[rel_window[0][1] - rel_window[0][0] + 1])
                                se_right = [lamb(element) for element in shanon_entropy[rel_window[0][1] -
                                                                                        rel_window[0][0] + 2:
                                                                                        max_window]]
                                print("%sshanon entropy%s:%s%s%s%s, %s, %s%s%s\n "
                                      % (white, " " * (space[0] - len("shanon entropy")),
                                         end, blue, str(se_left)[1:-1], end, str(se_phospho),
                                         green, str(se_right)[1:-1], end))
                                print("%sinformation content :%s%s%s%s%s%s,%s,%s%s%s%s\n "
                                      % (white, " " * (int(space[1] / 2) - 3), end, blue,
                                         lamb(IC[0]), " " * (int(space[1] / 2) - 3), end,
                                         lamb(IC[2]), " " * (int(space[1] / 2) - 3),
                                         green, lamb(IC[1]), end))
                                print("ACH alignment%s:%s%s%s%s%s,%s,%s%s%s%s\n "
                                      % (" " * (space[0] - len("ACH alignment")), blue,
                                         " " * (int(space[1] / 2) - 3), lamb(ACH[0]), end,
                                         " " * (int(space[1] / 2) - 3), lamb(ACH[2]), green,
                                         " " * (int(space[1] / 2) - 3), lamb(ACH[1]), end))
                                print("ACH prot%s:%s%s%s%s%s,%s,%s%s%s%s \n "
                                      % (" " * (space[0] - len("ACH prot")), blue,
                                         " " * (int(space[1] / 2) - 3), lamb(ACH_prot[0]), end,
                                         " " * (int(space[1] / 2) - 3), lamb(ACH_prot[2]), green,
                                         " " * (int(space[1] / 2) - 3), lamb(ACH_prot[1]), end))


def find_neg(phospho_sites, pattern, max_window, color, nthread, file, align_ortho_window):
    path = "%s/data" % os.path.abspath(os.path.dirname
                                       (os.path.dirname(__file__)))
    path2align = '%s/align' % path
    with open(file, 'a', newline='') as g:
        writer = csv.writer(g, delimiter=";")
        data_thread = []
        thread_list =[]
        for i in range(1, nthread + 1):
            data_thread.append({key: value for j, (key, value) in enumerate(phospho_sites.items())
                                if j % i == i-1})
        for data in data_thread:
            thread_list.append(fill_neg_sites(pattern, data, path2align, max_window, color, path,
                                              writer, align_ortho_window))
        for thread in thread_list:
            thread.start()
        for thread in thread_list:
            thread.join()

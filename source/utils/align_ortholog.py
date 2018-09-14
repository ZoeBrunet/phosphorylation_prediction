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
from Bio import SeqIO
from biothings_client import get_client
from difflib import SequenceMatcher
from source.utils.sequence import is_metazoan


def run_muscle(file_input):
    file = os.path.basename(file_input)
    file_output = "%s_align.fasta" % os.path.splitext(file)[0]
    path = os.path.dirname(os.path.dirname(file_input))
    path2align = "%s/align" % path
    if not os.path.exists(path2align):
        os.mkdir(path2align)
    path2outfile = "%s/%s" % (path2align, file_output)
    if not os.path.exists(path2outfile):
        muscle_cline = MuscleCommandline(input=file_input, out=path2outfile)
        os.system(str(muscle_cline))
    return(path2outfile)


def split_fasta(file):
    file_name = os.path.basename(file)
    path = os.path.dirname(os.path.dirname(file))
    os.makedirs("%s/metazoa" % path, exist_ok=True)
    os.makedirs("%s/non_metazoa" % path, exist_ok=True)
    path2metazoa = "%s/metazoa/%s_metazoa.fasta" % (path, file_name[:-6])
    path2nonmetazoa = "%s/non_metazoa/%s_non_metazoa.fasta" % (path, file_name[:-6])
    for record in SeqIO.parse(open(file), "fasta"):
        mt = get_client("taxon")
        position = str(record.id).find(":")
        taxID = record.id[:position]
        metazoa = is_metazoan(float(taxID), mt)
        f_out = path2metazoa if metazoa else path2nonmetazoa
        SeqIO.write([record], open(f_out, 'a'), "fasta")


def add_seq(record, f_out, taxid_list, taxid):
    r = SeqIO.write([record], open(f_out, 'a'), "fasta")
    if r != 1:
        print("Error while writing sequence: %s" % record.id)
    taxid_list.append(taxid)


def sort_fasta(path2cluster, taxID, sequence):
    taxid_list = []
    file = os.path.basename(path2cluster)
    path = os.path.dirname(os.path.dirname(path2cluster))
    os.makedirs("%s/sorted_fastas" % path, exist_ok=True)
    f_out = "%s/sorted_fastas/%s_sorted.fasta" % (path, file[:-6])
    seq_find = False
    if not os.path.exists(f_out):
        for record in SeqIO.parse(open(path2cluster), "fasta"):
            to_add = True
            taxonomy = str(record.id).split(":")[0]
            if float(taxonomy) == float(taxID) and taxID not in taxid_list and not seq_find:
                match = SequenceMatcher(None, record.seq,
                                        sequence).find_longest_match(0, len(record.seq),
                                                                     0, len(sequence))
                if match.size == 13:
                    add_seq(record, f_out, taxid_list, taxonomy)
                    seq_find = True
                else:
                    to_add = False
            if taxonomy not in taxid_list and to_add:
                add_seq(record, f_out, taxid_list, taxonomy)
    if not seq_find:
        os.remove(f_out)
        print("No corresponding sequence in fasta for %s" % file[:-6])
    else:
        path2align = run_muscle(f_out)
        split_fasta(path2align)

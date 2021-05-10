#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import re
import glob
from Bio import SeqIO

MSA_file = sorted([int(re.search("OCG(\d+)", file).group(1)) for file in glob.glob("*.cds.pal2nal.clean.pruned.long.muscle.trimal.fasta_12")])

OG_gene_seq_dict = {}
for OG in MSA_file:
    OG_gene_seq_dict[OG] = \
        {re.search("(.*)?\|", seq_record.id).group(1): str(seq_record.seq) for seq_record in
         SeqIO.parse("OCG%s.cds.pal2nal.clean.pruned.long.muscle.trimal.fasta_12" % OG, "fasta")}
#OG_gene_seq_dict = sorted(OG_gene_seq_dict)
#print OG_gene_seq_dict

# 获得物种名称，因为不是绝对单拷贝的数据，考虑到有的OGs会缺少个别物种，因此把所有OGs中的物种名都拿出来，然后去重，就能保证每个物种都不会漏
species_name_set = []
for seq_record in OG_gene_seq_dict.values():
    species_name_set += seq_record.keys()
species_name_set = sorted(set(species_name_set))
#print species_name_set, len(species_name_set)

#输出fasta格式的单拷贝orthogroup的supermatrix
with open("supermatrix_%s_OG_codon12.fasta" % len(MSA_file), "w") as fw:
    for sp in species_name_set:
        sequence = ""
        for og, seq_record in sorted(OG_gene_seq_dict.items()):
            length = [len(seq_record[i]) for i in species_name_set if i in seq_record.keys()][0]
#            print sp + "\t" + "OCG" + str(og) + "\t" + str(length)
            if sp in seq_record.keys():
                sequence += seq_record[sp]
            else:
                sequence += length * "-"
        fw.write(">" + sp + "\n" + sequence + "\n")

# 输出phlip格式的单拷贝orthogroup的supermatrix
with open("supermatrix_%s_OG_codon12.phy" % len(MSA_file), "w") as fw:
    # header
    length = 0
    for og, seq_record in sorted(OG_gene_seq_dict.items()):
        seq_len = len(seq_record.values()[0])
        length += seq_len
    fw.write(str(len(species_name_set)) + "\t" + str(length) + "\n")
    # sequence
    for sp in species_name_set:
        sequence = ""
        for og, seq_record in sorted(OG_gene_seq_dict.items()):
            length = [len(seq_record[species]) for species in species_name_set if species in seq_record.keys()][0]
            if sp in seq_record.keys():
                sequence += seq_record[sp]
            else:
                sequence += length * "-"
        fw.write(sp + "\t" + sequence + "\n")

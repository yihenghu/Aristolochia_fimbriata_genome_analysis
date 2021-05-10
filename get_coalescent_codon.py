import glob
import re
from Bio import SeqIO


def rm_3(seq):
    new = ''
    for index, i in enumerate(str(seq)):
        if (index + 1) % 3 != 0:
            new += i
    return new


def save_3(seq):
    new = ''
    for index, i in enumerate(str(seq)):
        if (index + 1) % 3 == 0:
            new += i
    return new


fileID_list = []
for file in glob.glob('*.fasta'):
    fileID = int(
        re.search(
            'OCG(\d+).cds.pal2nal.clean.pruned.long.muscle.trimal.fasta',
            file).group(1))
    fileID_list.append(fileID)

# Coalescent-based tree estimated from concatenated alignments of single gene first codon and second codon positions
for fileID in sorted(fileID_list):
    ID2seq = {}
    seqRs = SeqIO.parse('OCG%s.cds.pal2nal.clean.pruned.long.muscle.trimal.fasta' % fileID, 'fasta')
    for seqR in seqRs:
        ID2seq[seqR.id] = rm_3(seqR.seq)

    with open('OCG%s.cds.pal2nal.clean.pruned.long.muscle.trimal.fasta_12' % fileID, 'w') as fw:
        for ID, seq in ID2seq.items():
            fw.write('>' + ID + '\n' + seq + '\n')

# Coalescent-based tree estimated from concatenated alignments of single gene third codon positions
for fileID in sorted(fileID_list):
    ID2seq = {}
    seqRs = SeqIO.parse('OCG%s.cds.pal2nal.clean.pruned.long.muscle.trimal.fasta' % fileID, 'fasta')
    for seqR in seqRs:
        ID2seq[seqR.id] = save_3(seqR.seq)

    with open('OCG%s.cds.pal2nal.clean.pruned.long.muscle.trimal.fasta_3' % fileID, 'w') as fw:
        for ID, seq in ID2seq.items():
            fw.write('>' + ID + '\n' + seq + '\n')

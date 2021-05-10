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
for file in glob.glob('*fa'):
    fileID = int(
        re.search(
            'OCG(\d+).mafftHA.pruned.gb.filter.fa',
            file).group(1))
    fileID_list.append(fileID)


# concatenated alignments of first and second codon positions
sp2seq = {}

for fileID in sorted(fileID_list):
    seqRs = SeqIO.parse('OCG%s.mafftHA.pruned.gb.filter.fa' % fileID, 'fasta')
    for seqR in seqRs:
        sp = re.search('(\w*)|', seqR.id).group(1)
        if sp in sp2seq:
            sp2seq[sp] += rm_3(seqR.seq)
        else:
            sp2seq[sp] = rm_3(seqR.seq)

with open('all.fa.12', 'w') as fw:
    for sp, seq in sp2seq.items():
        fw.write('>' + sp + '\n' + seq + '\n')

# concatenated alignments of third codon positions
sp2seq = {}

for fileID in sorted(fileID_list):
    seqRs = SeqIO.parse('OCG%s.mafftHA.pruned.gb.filter.fa' % fileID, 'fasta')
    for seqR in seqRs:
        sp = re.search('(\w*)|', seqR.id).group(1)
        if sp in sp2seq:
            sp2seq[sp] += save_3(seqR.seq)
        else:
            sp2seq[sp] = save_3(seqR.seq)

with open('all.fa.3', 'w') as fw:
    for sp, seq in sp2seq.items():
        fw.write('>' + sp + '\n' + seq + '\n')

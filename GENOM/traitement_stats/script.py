from Bio import SeqIO
from Bio import AlignIO
import os
from Bio.Align.Applications import ClustalOmegaCommandline
from itertools import combinations
from Bio.SeqUtils import GC, GC_skew
from Bio import pairwise2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import re
def list_dir(basepath):
    fastas = []
    for entry in os.listdir(basepath):
        fastas.append(os.path.join(basepath, entry))
    return fastas


def read_sequences(filepath):
    sequences =[]
    for seq_record in SeqIO.parse(filepath, "fasta"):
        sequences.append(seq_record)
    return sequences



def construct_alignment(fastas):
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))


def write_sequences(filepath,sequences):
    with open(filepath, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')


def align_sequences(infile,outfile):
    #alignments = pairwise2.align.globalxx(seq1, seq2)
    command = './files/clustalo/clustalo -i '+infile +' -o '+ outfile +' --outfmt=phy --force'
    os.system(command)

# def read_alignment(filename):
#     align = AlignIO.read(filename, "phylip")
#     return align
#
# def write_alignement(filepath,alignement):
#     with open(filepath, 'w') as handle:
#         AlignIO.write(alignement, handle, 'phylip')

def select_sequence(sequences, species, ensure_all=False):
    selected_sequences =[]
    for sp in species:
        for seq in sequences:
            if sp in seq.id:
                seq.id = sp
                selected_sequences.append(seq)
                break
    if len(selected_sequences) != len(species) and ensure_all:
        return None
    return selected_sequences

def read_gff(filepath):
    file = open(filepath,"r")
    for line in file.readlines():
        line=line.split("\t")
        debut = line[3]
        fin = line[4]
        regex = r"Gene=(.+?);"
        name = re.search(regex, line[-1]).group(1)
    file.close()

if __name__ == '__main__':
    input_path = "../sequences_geniques/"
    trg = {}
    bayanus = read_sequences(input_path + "Sbayanus.fsa" )
    cerevisiae = read_sequences(input_path + "Scerevisiae.fsa" )
    kudriavzevii = read_sequences(input_path + "Skudriavzevii.fsa" )
    mikatae = read_sequences(input_path + "Smikatae.fsa" )
    paradoxus = read_sequences(input_path + "Sparadoxus.fsa" )
    read_gff(input_path + "Sbayanus.gff")

    # gc_stat = []
    # length_stat = []
    # # gc_skew_stat = []
    # for record in bayanus:
    #     gc_stat.append(GC(record.seq))
    #     if len(record) > 10000:
    #         print(len(record))
    #         print(record)
    #     length_stat.append(np.log(len(record)))
    #     # gc_skew_stat.append(GC_skew(record.seq))
    # plt.figure()
    # plt.hist(gc_stat,100)
    # plt.figure()
    # plt.hist(length_stat,100, weights=np.ones(len(length_stat)) / len(length_stat))
    # plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    # plt.show()

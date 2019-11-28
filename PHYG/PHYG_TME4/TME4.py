from Bio import SeqIO
import os
from Bio.Align.Applications import ClustalOmegaCommandline
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Exercise 2 :                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

def list_dir(basepath):

    fastas = []
    for entry in os.listdir(basepath):
        if os.path.isfile(os.path.join(basepath, entry)):
            fastas.append(os.path.join(basepath, entry))
    return fastas

def construct_alignment(fastas):
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

def extract_sequences(filepath):
    return SeqIO.parse(filepath, "fasta")

def align_sequences(sequences):
    clustalomega = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)


if __name__ == '__main__':
    print(extract_sequences('files/TME4_sequences/PF01599.fasta'))
    with open("files/species.list", 'r') as f:
        species = f.read().splitlines()


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Exercise 3 :                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Exercise 4 :                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

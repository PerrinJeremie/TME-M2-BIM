from Bio import SeqIO
from Bio import AlignIO
import os
from Bio.Data import CodonTable
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align.Applications import ClustalOmegaCommandline
from itertools import combinations
from Bio.SeqUtils import GC, GC_skew
from Bio import pairwise2
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import codonalign
from Bio.Seq import translate
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.codonalign import build
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from scipy import stats
import re
import itertools

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


def write_sequences(filepath,sequences):
    with open(filepath, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')

def read_trg(filepath):
    file = open(filepath,"r")
    txt = file.read().split("\n")
    trg=[]
    for line in txt:
        trg.append(line.split("\t"))
    file.close()
    return trg

def read_gff(filepath):
    positions = {}
    file = open(filepath,"r")
    for line in file.readlines():
        line=line.split("\t")
        start = int(line[3])
        end = int(line[4])
        regex = r"Gene=(.+?);"
        name = re.search(regex, line[-1]).group(1)
        positions[name] = {"start":start,"end":end}
    file.close()
    return positions

def test_trg(gene,trgs,trg_seq):
    for trg in trgs:
        if str(gene.id) in trg:
            trg_seq[str(gene.id)] = gene
            return True
    return False
def test_orfan(orphans,gene):
    return str(gene.id) in orphans

def distance_to_nearest(gene,position):
    prefix,num = gene.split('.')
    num = int(num)
    gene_start = position[gene]["start"]
    gene_end = position[gene]["end"]
    upstream_dist = float('inf')
    k = num
    while k > 1:
        k -= 1
        if prefix + "." + str(k) in positions:
            upstream_dist = gene_start - positions[prefix + "." + str(k)]["end"]
            break
    k = num
    downstream_dist = float('inf')
    while k < 5000:
        k += 1
        if prefix + "." + str(k) in positions:
            downstream_dist = positions[prefix + "." + str(k)]["start"] - gene_end
            break
    distance = min(upstream_dist,downstream_dist)

    if distance == float('inf'):
        return -1
    return distance
def compute_identity(gene1,gene2):
    aln = pairwise2.align.globalxx(gene1 ,gene2)[0]
    matches = sum(aa1 == aa2 for aa1, aa2 in zip(aln[0], aln[1]))
    pct_identity = 100.0 * matches / len(aln[0])
    return pct_identity
def compute_kaks(gene1,gene2):

    prot1 = str(translate(gene1.seq))
    prot2 = str(translate(gene2.seq))

    # matrix = matlist.blosum62
    # gap_open = -10
    # gap_extend = -0.5
    # aln = pairwise2.align.globalds(prot1 ,prot2, matrix, gap_open, gap_extend)[0]
    aln = pairwise2.align.globalxx(prot1 ,prot2)[0]
    prot1, prot2, score, begin, end = aln

    nuc1 = SeqRecord(Seq(str(gene1.seq), alphabet=IUPAC.IUPACUnambiguousDNA()), id='nuc1')
    nuc2 = SeqRecord(Seq(str(gene2.seq), alphabet=IUPAC.IUPACUnambiguousDNA()), id='nuc2')
    prot1 = SeqRecord(Seq(prot1, alphabet=IUPAC.protein),id='pro1')
    prot2 = SeqRecord(Seq(prot2, alphabet=IUPAC.protein),id='pro2')

    aln = MultipleSeqAlignment([prot1, prot2])

    codon_aln = codonalign.build(aln, [nuc1, nuc2],max_score=10)

    dn,ds = cal_dn_ds(codon_aln[0], codon_aln[1], method="ML",codon_table=CodonTable.unambiguous_dna_by_id[1])
    print(dn,ds,len(gene1))
    return dn/ds

if __name__ == '__main__':
    input_path = "../sequences_geniques/"
    trg_seq = {}
    orphans = np.genfromtxt("orphan_final2",dtype=str)
    trgs = read_trg("trgs_final2")
    gc_stat = []
    length_stat = []
    distance_stat = []

    gc_stat_orf = []
    length_stat_orf = []
    distance_stat_orf = []

    gc_stat_trg = []
    length_stat_trg = []
    distance_stat_trg = []


    for species in {"Sbayanus","Scerevisiae","Skudriavzevii","Smikatae","Sparadoxus"}:
        genome = read_sequences(input_path + species + ".fsa" )
        positions = read_gff(input_path + species + ".gff")
        for record in genome:
            distance = distance_to_nearest(str(record.id),positions)
            if test_trg(record,trgs,trg_seq):
                gc_stat_trg.append(GC(record.seq))
                length_stat_trg.append(np.log(len(record))/np.log(2))
                if distance >= 0:
                    if distance == 0:
                        distance_stat_trg.append(0)
                    else:
                        distance_stat_trg.append(np.log(distance)/np.log(10))
            elif test_orfan(orphans,record):
                gc_stat_orf.append(GC(record.seq))
                length_stat_orf.append(np.log(len(record))/np.log(2))
                if distance >= 0:
                    if distance == 0:
                        distance_stat_orf.append(0)
                    else:
                        distance_stat_orf.append(np.log(distance)/np.log(10))
            else:
                if int(re.search(r"LEN:(.+?);",record.description).group(1)) > 0:
                    gc_stat.append(GC(record.seq))
                    length_stat.append(np.log(len(record))/np.log(2))
                    if distance >= 0:
                        if distance == 0:
                            distance_stat.append(0)
                        else:
                            distance_stat.append(np.log(distance)/np.log(10))

    # identity_stat = []
    # for trg in trgs:
    #     identity_temp = []
    #     for pair in itertools.combinations(trg,2):
    #         identity_temp = compute_identity(trg_seq[pair[0]],trg_seq[pair[1]])
    #     identity_stat.append(np.mean(np.array(identity_temp)))

    # print(compute_kaks(trg_seq[trg[0][0]],trg_seq[trg[0][1]]))

    plt.figure()
    plt.hist(gc_stat,100,weights=np.ones(len(gc_stat)) / len(gc_stat))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("GC content (%)")
    plt.title("GC content of non de novo gene")
    plt.suptitle("Nb gene = " + str(len(gc_stat)),fontsize=8)

    plt.figure()
    plt.hist(gc_stat_orf,100,weights=np.ones(len(gc_stat_orf)) / len(gc_stat_orf))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("GC content (%)")
    plt.title("GC content of orphan gene")
    plt.suptitle("Nb gene = " + str(len(gc_stat_orf)),fontsize=8)

    plt.figure()
    plt.hist(gc_stat_trg,100,weights=np.ones(len(gc_stat_trg)) / len(gc_stat_trg))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("GC content (%)")
    plt.title("GC content of TRG gene")
    plt.suptitle("Nb gene = " + str(len(gc_stat_trg)),fontsize=8)

    print("CG content of gene")
    print("Classic genes and orphans")
    t, p = stats.ttest_ind(gc_stat,gc_stat_orf,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))
    print("Classic genes and TRG")
    t, p = stats.ttest_ind(gc_stat,gc_stat_trg,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))
    print("TRG genes and orphans")
    t, p = stats.ttest_ind(gc_stat_trg,gc_stat_orf,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))

    plt.show()

    plt.figure()
    plt.hist(length_stat,100, weights=np.ones(len(length_stat)) / len(length_stat))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("log2(size)")
    plt.title("Gene length of non de novo gene")
    plt.suptitle("Nb gene = " + str(len(length_stat)),fontsize=8)

    plt.figure()
    plt.hist(length_stat_orf,100, weights=np.ones(len(length_stat_orf)) / len(length_stat_orf))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("log2(size)")
    plt.title("Gene length of orphan gene")
    plt.suptitle("Nb gene = " + str(len(length_stat_orf)),fontsize=8)

    plt.figure()
    plt.hist(length_stat_trg,100, weights=np.ones(len(length_stat_trg)) / len(length_stat_trg))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("log2(size)")
    plt.title("Gene length of TRG gene")
    plt.suptitle("Nb gene = " + str(len(length_stat_trg)),fontsize=8)

    print("Length of gene")
    print("Classic genes and orphans")
    t, p = stats.ttest_ind(length_stat,length_stat_orf,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))
    print("Classic genes and TRG")
    t, p = stats.ttest_ind(length_stat,length_stat_trg,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))
    print("TRG genes and orphans")
    t, p = stats.ttest_ind(length_stat_trg,length_stat_orf,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))

    plt.show()

    plt.figure()
    plt.hist(distance_stat,100, weights=np.ones(len(distance_stat)) / len(distance_stat))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("log10(distance to nearest gene)")
    plt.title("Distance to nearest gene for non de novo gene")
    plt.suptitle("Nb gene = " + str(len(distance_stat)),fontsize=8)

    plt.figure()
    plt.hist(distance_stat_orf,100, weights=np.ones(len(distance_stat_orf)) / len(distance_stat_orf))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("log10(distance to nearest gene)")
    plt.title("Distance to nearest gene for orphan gene")
    plt.suptitle("Nb gene = " + str(len(distance_stat_orf)),fontsize=8)

    plt.figure()
    plt.hist(distance_stat_trg,100, weights=np.ones(len(distance_stat_trg)) / len(distance_stat_trg))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel("log10(distance to nearest gene)")
    plt.title("Distance to nearest gene for TRG gene")
    plt.suptitle("Nb gene = " + str(len(distance_stat_trg)),fontsize=8)



    print("Distance to nearest gene")
    print("Classic genes and orphans")
    t, p = stats.ttest_ind(distance_stat,distance_stat_orf,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))
    print("Classic genes and TRG")
    t, p = stats.ttest_ind(distance_stat,distance_stat_trg,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))
    print("TRG genes and orphans")
    t, p = stats.ttest_ind(distance_stat_trg,distance_stat_orf,equal_var=False)
    print("t = " + str(t))
    print("p = " + str(p))

    # plt.figure()
    # plt.hist(identity_stat,100,weights=np.ones(len(identity_stat)) / len(identity_stat))
    # plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    # plt.xlabel("Identity (%)")
    # plt.title("Identity between trg groups")
    # plt.suptitle("Nb trg group = " + str(len(identity_stat)),fontsize=8)
    plt.show()

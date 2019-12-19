from Bio import SeqIO
from Bio import AlignIO
import os
from Bio.Align.Applications import ClustalOmegaCommandline
from itertools import combinations
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Exercise 2 :                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

def construct_alignment(fastas):
    for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

def read_sequences(filepath):
    sequences =[]
    for seq_record in SeqIO.parse(filepath, "fasta"):
        sequences.append(seq_record)
    return sequences

def write_sequences(filepath,sequences):
    with open(filepath, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')


def align_sequences(infile,outfile):
    command = './files/clustalo/clustalo -i '+infile +' -o '+ outfile +' --outfmt=phy --force'
    os.system(command)

def read_alignment(filename):
    align = AlignIO.read(filename, "phylip")
    return align

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

def create_distance_matrix(infile,outfile):
    command = "fprotdist " + infile +" "+ outfile
    os.system(command)

def build_tree_NJ(infile,outfile):
    command = "fneighbor " + infile +" outputs/output -outtreefile "  + outfile
    os.system(command)
    os.system("rm outputs/output")
    os.system("rm "+infile)

def build_tree_ML(infile,outfile):
    command = "fproml " + infile +" -outfile outputs/output -intreefile '' -outtreefile "  + outfile
    os.system(command)
    os.system("rm outputs/output")

# if __name__ == '__main__':
#     output_path = "outputs/ex2/"
#     sequences = read_sequences('files/TME4_sequences/PF01599.fasta')
#     with open("files/species.list", 'r') as f:
#         species = f.read().splitlines()
#     selected_sequences = select_sequence(sequences, species)
#     write_sequences(output_path +'species_ex2.fasta', selected_sequences)
#     align_sequences(output_path +'species_ex2.fasta',output_path +'alignment_ex2.phylip')
#     create_distance_matrix(output_path +'alignment_ex2.phylip',output_path +"mat.txt")
#     build_tree_NJ(output_path +"mat.txt",output_path + "NJ_ex2.tree")
#     build_tree_ML(output_path +'alignment_ex2.phylip',output_path +"ML_ex2.tree")




##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Exercise 3 :                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

def list_dir(basepath):
    fastas = []
    for entry in os.listdir(basepath):
        fastas.append(os.path.join(basepath, entry))
    return fastas

def write_alignement(filepath,alignement):
    with open(filepath, 'w') as handle:
        AlignIO.write(alignement, handle, 'phylip')



# if __name__ == '__main__':
#     basepath = 'files/TME4_sequences/'
#     output_path = "outputs/ex3/"
#     fastas = list_dir(basepath)
#     with open("files/species.list", 'r') as f:
#         species = f.read().splitlines()
#     i=0
#     for fa in fastas:
#         sequences = read_sequences(fa)
#         selected_sequences =select_sequence(sequences, species,True)
#         if selected_sequences == None:
#             continue
#         write_sequences(output_path +'temp_species_ex3.fasta', selected_sequences)
#         align_sequences(output_path +'temp_species_ex3.fasta',output_path +'temp_alignment_ex3.phylip')
#         alignment = read_alignment(output_path +'temp_alignment_ex3.phylip')
#         alignment.sort()
#         if i == 0:
#             concatenate_alignement = alignment
#             i+=1
#         else:
#             concatenate_alignement += alignment
#     write_alignement(output_path +"concatenate_alignement_exo3.phylip",concatenate_alignement)
#     create_distance_matrix(output_path +'concatenate_alignement_exo3.phylip',output_path +"mat.txt")
#     build_tree_NJ(output_path +"mat.txt",output_path + "NJ_ex3.tree")
#     build_tree_ML(output_path +'concatenate_alignement_exo3.phylip',output_path +"ML_ex3.tree")
#     os.system("rm outputs/ex3/*temp*")



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Exercise 4 :                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

def read_clades(infile):
    species_by_clade = {}
    with open(infile, 'r') as f:
             lines = f.read().splitlines()
    for line in lines:
        if '#' in line:
            clade = line.replace('#',"").lower()
            species_by_clade[clade] = []
        else:
            species_by_clade[clade].append(line)
    return species_by_clade


if __name__ == '__main__':
    output_path = "outputs/ex4/"
    data = {}
    species_by_clade = read_clades("files/clades.list")
    basepath = 'files/TME4_clades'
    famillies = list_dir(basepath)
    for familly in famillies:
        familly_name = familly.split("/")[-1]
        fastas = list_dir(familly)
        all = True
        temp = {}
        for fa in fastas:
            clade_name = fa.split(".")[-2].split('/')[-1]
            if clade_name not in species_by_clade:
                continue
            sequences = read_sequences(fa)
            selected_sequences = select_sequence(sequences, species_by_clade[clade_name], ensure_all=True)
            if selected_sequences == None :
                all = False
                break
            temp[clade_name] = selected_sequences
        if all:
            data[familly_name] = temp

    clades_pairs = combinations(species_by_clade.keys(),2)
    for pair in clades_pairs:

        first = True
        clade1 = pair[0]
        clade2 = pair[1]
        output =  output_path + clade1 + "-" + clade2
        for familly in data.keys():

            try:
                write_sequences(output_path + 'temp.fasta', data[familly][clade1] + data[familly][clade2])
            except:
                continue
            align_sequences(output_path+'temp.fasta',output_path+'temp.phylip')
            alignment = read_alignment(output_path+'temp.phylip')
            alignment.sort()
            if first:
                concatenate_alignement = alignment
                first = False
            else:
                concatenate_alignement += alignment
        write_alignement(output+ ".phylip",concatenate_alignement)
        create_distance_matrix(output+ ".phylip",output_path +"mat.txt")
        build_tree_NJ(output_path + "mat.txt",output + "_NJ.tree")
        build_tree_ML(output + ".phylip" , output +"_ML.tree")
    os.system("rm outputs/ex4/*temp*")

# -*- coding: Utf-8 -*-

# developped by Aubin Fleiss
# contact : aubin.fleiss@gmail.com

import random as rd
import numpy as np
import copy
import sys

""" This function return a non-zero value according to a Poisson distribution
with the lambda parameter provided in input """
def poisson(lam):
	val = 0
	while val == 0:
		val = int(np.random.poisson(lam, 1))
	return val

""" This function takes in input a genome and returns a pair where
the first value is a chromosome identifier and the second value is a random
index within the list of genes of the choosen chromosome """
def choose_coordinates(genome):
	# choose chromosome : adjust probability based on chromosome length
	chroms=genome.keys()
	nbGenes=sum([len(genome[chrom]) for chrom in genome])
	p=[float(len(genome[chrom]))/float(nbGenes) for chrom in sorted(genome.keys())]
	chromChoice = int(np.random.choice(sorted(genome.keys()), 1, p=p))
	# choose gene in previously selected chromosome
	geneChoice = rd.choice(range(len(genome[chromChoice])))
	return((chromChoice,geneChoice))

""" This function takes in input a genome and prints it to the standard output """
def print_genome(genome):
	for chrom in genome.keys():
		print(chrom),
		print(genome[chrom])

""" EXERCISE 2 - Inversion function
 	This function takes an input genome and return it after applying an inversion event
	It use a mean of the length to determine the number of genes that will be inverted"""
def inversion(genome,mean_inv_len):
	chrom_index, gene_index = choose_coordinates(genome)
	inversion_len = poisson(mean_inv_len)
	while len(genome[chrom_index]) <= gene_index + inversion_len :
		chrom_index, gene_index = choose_coordinates(genome)
	genome[chrom_index][gene_index:gene_index+inversion_len] = [-x for x in genome[chrom_index][gene_index:gene_index+inversion_len][::-1]]
	return(genome)

""" EXERCISE 2 - Deletion function
	This function takes an input genome and return it after applying an deletion event
	It use a mean of deletion length to determine the number of genes that will be deleted"""
def deletion(genome,mean_del_len):
	chrom_index, gene_index = choose_coordinates(genome)
	delete_len = poisson(mean_del_len)
	while delete_len > len(genome[chrom_index])-gene_index:
		delete_len = poisson(mean_del_len)
	if 0 in genome[chrom_index][gene_index:gene_index+delete_len]:
		genome[chrom_index][gene_index:gene_index+delete_len+1] = [0]
	else:
		genome[chrom_index][gene_index:gene_index+delete_len] = []
	return(genome)

""" EXERCISE 2 - Fission function
	This function takes an input genome and return it after applying an fission event"""
def fission(genome):
	chrom_index, gene_index = choose_coordinates(genome)
	while gene_index == 0 :
		chrom_index, gene_index = choose_coordinates(genome)
	genome[max(genome.keys()) + 1] = genome[chrom_index][gene_index:]
	genome[chrom_index] = genome[chrom_index][:gene_index]
	if 0 in genome[chrom_index]:
		genome[max(genome.keys())].insert(rd.randint(0,len(genome[max(genome.keys())])),0)
	else:
		genome[chrom_index].insert(rd.randint(0,len(genome[chrom_index])),0)

	if (not 0 in genome[chrom_index]) or (not 0 in genome[max(genome.keys())]):
		print(genome[chrom_index])
		print(genome[max(genome.keys())]) 
	return(genome)


def main( argv=None ):

	# definition of a test genome
	genome = {}
	genome[1]=range(1,22)
	genome[1].insert(5,0)
	genome[2]=range(22,54)
	genome[2].insert(10,0)
	genome[3]=range(54,76)
	genome[3].insert(7,0)
	genome[4]=range(76,101)
	genome[4].insert(11,0)

	print("### ancestor genome ###")
	print_genome(genome)

	print("### genome after an inversion ###")
	genome = inversion(genome,5)
	print_genome(genome)

	print("### genome after a deletion ###")
	genome=deletion(genome,1)
	print_genome(genome)

	print("### genome fission ###")
	genome=fission(genome)
	print_genome(genome)

	return 0


if __name__=="__main__":
	sys.exit(main())

 #											-*- Mode: Python -*-
#											-*- coding UTF-8 -*-
# topsequenceoverlap.py
# Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
# Author :  Adam-Nicolas Pelletier
# Last modified On: 2015 - 04- 13

from __future__ import generators
import numpy as np
import pandas as pd
import itertools
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import math
import os
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()



########################################################################################################################################################
########################################################## USER INPUT ##################################################################################
## GOAL: Generate a synthetic promoter based on a collection of unified PWMs (pwmunif.py) from their top sequence. 
## The top sequence is defined as the sequence with the higest probability from the PWM

seq = "ATGCTGA"
file1 = "Homosapiens_CISBP/pwmsLOG/M0186_1.01.txt"
file2 = "Homosapiens_CISBP/pwmsLOG/M0374_1.01.txt"

pgl4RE = [NheI, XhoI, HindIII] #list of RE in order on PGL4 vector, for reference purposes. EcoRV was removed: not compatible with Renilla. 
pgl4REstr = []
for i in list(pgl4RE):
	pgl4REstr.append(str(i))

fivepr = "AATTAT"  # write the sequence you want to add in the 5' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
threepr = "TGTAAC" # write the sequence you want to add in the 3' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
atrich = "ATATCAGAAT"


longpwm = 24
########################################################################################################################################################
########################################################################################################################################################

dic = {0:"A",1:"C",2:"G",3:"T"}  # dictionary of indexes of each nucleotide for matrices


def randomseq(length):   # generate a random n length DNA sequence as a matrix
	matrix = np.zeros( (4, length) )
	index = []
	for i in range (length):    
	    index.append([random.randrange(0,4,1), i]) 
	a = np.array(index)  
	matrix[a[:,0], a[:,1]] = 1
	return matrix

def matrixconverter(seqmatrix, dic): #generates a DNA sequence string based on a NumPy matrix.
	a = np.transpose(np.nonzero(np.transpose(seqmatrix))).tolist()
	seqstring = ""
	for i in a:
		seqstring += dic[i[1]]
	return seqstring

def topseq(filename):  
	tf1 = np.loadtxt(filename, skiprows=1)
	tf2 = tf1[:,1:].transpose()
	indexes = np.argmax(tf2,axis=0).tolist()  #generates a list of indexes of the highest probability base per position. 
	seqstring = ""
	for i in indexes: #generates a string nucleotide sequence corresponding top nucleotide for each position
		seqstring += dic[i]
	return seqstring # returns the top motif


def remove_duplicates(values):
	output = []
	seen = set()
	for value in values:
		# If value has not been encountered yet,
		# ... add it to both list and set.
		if value not in seen:
			output.append(value)
			seen.add(value)
	return output




# Knuth-Morris-Pratt string matching
# David Eppstein, UC Irvine, 1 Mar 2002
#from http://code.activestate.com/recipes/117214/
def KnuthMorrisPratt(text, pattern):
 
    '''Yields all starting positions of copies of the pattern in the text.
Calling conventions are similar to string.find, but its arguments can be
lists or iterators, not just strings, it returns all matches, not just
the first one, and it does not need the whole text in memory at once.
Whenever it yields, it will have read the text exactly up to and including
the match that caused the yield.'''
 
    # allow indexing into pattern and protect against change during yield
    pattern = list(pattern)
 
    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift
 
    # do the actual search
    startPos = 0
    matchLen = 0
    for c in text:
        while matchLen == len(pattern) or \
              matchLen >= 0 and pattern[matchLen] != c:
            startPos += shifts[matchLen]
            matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == len(pattern):
            yield startPos

subsetlist = ["GATA2", "SPI1", "KLF3", "FOXO4", "HOXA6", "EGR2", "NFYA"]

filelist = []
for i in subsetlist:
	filelist.append("Homosapiens_CISBP/pwmsUnif"+"/"+i+".txt")

#### 1. Generate the list of PWM files to use
# filelist = []
# for i in os.listdir("Homosapiens_CISBP/pwmsUnif"):
# 	if i.endswith(".txt"):
# 		filelist.append("Homosapiens_CISBP/pwmsUnif"+"/"+i)

seqlist = []
for i in filelist:
	seqlist.append(topseq(i))





# Iterate possible versions of the promoter, and see which is the shortest. Loop rejects the sequences that are longeer than one previously generated. 

loop = 0
omegaseq = ""
omegalen = 2000


while loop < 5000:
	seqlistclean = remove_duplicates(seqlist)  # remove duplicate motifs
	fseq = " "
	resites = []
	random.shuffle(seqlistclean)		# here, the goal is to randomize the order of the mtifs in the final sequence. This allows to redo a sequence in the event that the first sequence generated 
	
									## is cleaved by our desired RE. Optional step. 
	fseq += seqlistclean[0]
	
	seqlistclean.pop(0)
	
	

	##Overlap information
	
	

	while len(seqlistclean) > 0:
		alphapos = 0 # position fo the overlap with the iterating sequence
		alphasize = -1 # length of the overlap of motif with iterating sequence. Start with -1 to allow for overlaps = 0 to be included towards the end of the loop. 
		alphaseq = ""  #seq of motif
		alphaindex = 0  #index of motif in motif lsit

		for i in seqlistclean:
			sequence = i
			d = list(KnuthMorrisPratt(fseq, sequence))

			if not d:
				betapos = len(sequence)
				betasize = 0
				betaindex = seqlistclean.index(sequence)
				betaseq = sequence
			else: 
				betapos = d[0]
				if len(sequence) + betapos > len(fseq):
					betasize = (len(sequence) + betapos) - len(fseq)
				else: 
					betasize = len(sequence)

				betaseq = sequence
				betaindex = seqlistclean.index(sequence)
				

			if betasize > alphasize:
				alphasize = betasize
				alphapos = betapos
				alphaseq = betaseq
				alphaindex = betaindex
				
		if alphasize == 0:
			fseq += alphaseq 
			seqlistclean.pop(alphaindex)

		elif (alphapos + len(alphaseq)) < len(fseq) - 1:
			seqlistclean.pop(alphaindex)
		else:
			fseq = fseq[0:len(fseq)] + alphaseq[alphasize:len(alphaseq)]
			seqlistclean.pop(alphaindex)


	if len(fseq) < omegalen:
		omegaseq = fseq
		omegalen = len(fseq)

	loop += 1  # this step was only added in the case where no RE combo could be used. It turns out most enzymes are very compatible with the set of motifs used. 

	print "Step: " + str(loop) + "   Length: " + str(omegalen)







resites= []

fseqf = Seq(omegaseq, amb)
resites.append(len(NheI.search(fseqf)))
resites.append(len(XhoI.search(fseqf)))
resites.append(len(HindIII.search(fseqf)))

dic = dict(zip(list(pgl4REstr), resites))
pgl4dict = dict(zip(list(pgl4REstr), list(pgl4RE)))

seqfinal = omegaseq.replace(" ","")



### Generate file


headerline = "Enzyme_pair" + "\t" + "Sequence" + "\n" 

with open("topsubset.txt", "w") as header:
	header.write(headerline)	

for i in list(itertools.combinations(list(pgl4REstr),2)):
	somme = [dic[list(i)[0]], dic[list(i)[1]]]
	if sum(somme) == 0:
		seqandsites = str(list(i)[0]) + "\t" + str(list(i)[1]) + "\n" + fivepr + pgl4dict[list(i)[0]].site + atrich + seqfinal + atrich + pgl4dict[list(i)[1]].site + threepr + "\n"
		output = str(seqandsites)
		with open("topsubset.txt", "a") as topsequence:
			topsequence.write(output)
	else: 
		pass

#### to add specific motifs at specific positions, we need to adapt the Knuth Morris approach. Snce all those specific positions are relative to the end, 
#and that we dont know the total length of the construct until the algorithm, we need to start writing it from the end. 
#TO achieve this, we need to revert (str[::-1]) every top sequence, and build it with those. Then, we integrate the motifs at specific positions FROM THE BEGGINING.
# Then we revert the whole sequence back once its complete. 
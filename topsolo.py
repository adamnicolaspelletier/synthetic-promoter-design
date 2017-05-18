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
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
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

pgl4RE = [KpnI, NheI, XhoI, HindIII] #list of RE in order on PGL4 vector, for reference purposes. EcoRV was removed: not compatible with Renilla. 
pgl4REstr = []
for i in list(pgl4RE):
	pgl4REstr.append(str(i))

fivepr = "AATTAT"  # write the sequence you want to add in the 5' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
threepr = "TGTAAC" # write the sequence you want to add in the 3' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
atrich = ["ATATCAGAAT","ATTCCATAGA"]
# atrich2 = "ATTCCATAGA"

###tfbs replicates
tfbsreplicates = 3


sp1 = "GGGGCGGGGC"
minimalp= "TAGAGGGTATATAATGGAAGCTCGACTTCCAG"

overlap = 25
longpwm = 24
########################################################################################################################################################
########################################################################################################################################################

dic = {0:"A",1:"C",2:"G",3:"T"}  # dictionary of indexes of each nucleotide for matrices

def reversecomp(rprimsequence): ## make a complement version of the sequence, and reverse it so it has the proper orientation
	a = ""
	tempzrev = rprimsequence
	tempzrev = tempzrev.replace("T","X")
	tempzrev = tempzrev.replace("A","T")
	tempzrev = tempzrev.replace("X","A")
	tempzrev = tempzrev.replace("C","Y")
	tempzrev = tempzrev.replace("G","C")
	tempzrev = tempzrev.replace("Y","G")
	templist = list(tempzrev)
	templist.reverse()
	for i in templist:
		a += i
	return a

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


coremotifdf = pd.read_csv("TF_core_update.txt", sep="\t", index_col = "Gene") 
core = coremotifdf[coremotifdf['Evidence'].isin(["D","I"])]

coredict = core["Motif"].to_dict()  # dict of Gne names as keys, core motif as values

filelist = []
seqdict = {}

del coredict["SP1"]
# subsetlist = ["GATA2", "SPI1", "KLF3", "FOXO4", "HOXA6", "EGR2", "NFYA"]

subsetlist = []
for i in coredict:
	subsetlist.append(i)

for i in xrange(len(subsetlist)):
	j = subsetlist[i]
	try: 
		filelist.append("Homosapiens_CISBP/pwmsUnif"+"/"+j+".txt")
		seqdict[j] = topseq("Homosapiens_CISBP/pwmsUnif"+"/"+j+".txt")
	except IOError:
		seqdict[j] = coredict[j]
		print 'IOError in file Homosapiens_CISBP/pwmsUnif/'+j+'.txt'

headerline = "Gene_Name" + "\t" + "Enzyme_1" + "\t" + "Enzyme_2" + "\t" + "Sequence" + "\t" + "Oligo_F" + "\t" + "Oligo_R" + "\n" 


outputname = "solo_promoters/oligosTFs_%s-replicates.txt" % str(tfbsreplicates)
with open(outputname, "w") as header:
		header.write(headerline)	
minimalseq = fivepr + pgl4RE[1].site + atrich[0] + minimalp + atrich[1] + pgl4RE[2].site + threepr
oligofmin = minimalseq[:((len(minimalseq)/2)+(overlap/2))]
oligormin = reversecomp(minimalseq[((len(minimalseq)/2)-(overlap/2)):])

mini = "Minimal_prom" + "\t"+ "NheI" + "\t"+ "XhoI" + "\t" + str(minimalseq) + "\t"+ str(oligofmin) + "\t"+ str(oligormin) + "\n"
with open(outputname, "a") as header:
		header.write(mini)

print seqdict

atlisttemp = tfbsreplicates*atrich
for j in seqdict:
	lencheck = 0
	
	atlist = atlisttemp


	while lencheck == 0:

		# outputname= "solo_promoters/%s.txt" % j
		# with open(outputname, "w") as header:
		# 	header.write(headerline)	

		resites= []
		omegaseq = seqdict[j]
		fseqf = Seq(omegaseq, amb)
		resites.append(len(KpnI.search(fseqf)))
		resites.append(len(NheI.search(fseqf)))
		resites.append(len(XhoI.search(fseqf)))
		resites.append(len(HindIII.search(fseqf)))
	  #NEED TO EXPORT ONLY THE Kpn1 / NheI.  XHOI sequence, it works with everything. and then design a Nhe1 xhoI minimal promoter site
		dic = dict(zip(list(pgl4REstr), resites))
		pgl4dict = dict(zip(list(pgl4REstr), list(pgl4RE)))
		

		
		seqtemp = omegaseq.replace(" ","")
		seqfinal = atlist[1]

		for i in xrange(tfbsreplicates):
			seqfinal += (seqtemp + atlist[i])
		

		

		seqandsites = fivepr + pgl4RE[0].site + seqfinal + pgl4RE[2].site + pgl4RE[1].site + threepr
		seq = str(seqandsites)



		lenseq = len(seqandsites)
		oligof = seqandsites[:((lenseq/2)+(overlap/2))]
		oligor = reversecomp(seqandsites[((lenseq/2)-(overlap/2)):])

		if len(oligof) <= 65 or len(oligor) <= 65:
			output = str(j) + "\t" + pgl4REstr[0] + "\t" + pgl4REstr[1] + "\t" + str(seqandsites) + "\t" + str(oligof) + "\t" + str(oligor) +"\n"
			
			with open(outputname, "a") as topsequence:
				topsequence.write(output)

			lencheck = 1

		elif len(atlist[0]) <= 5:
		
			output = str(j) + "\t" + pgl4REstr[0] + "\t" + pgl4REstr[1] + "\t" + str(seqandsites) + "\t" + str(oligof) + "\t" + str(oligor) +"\n"
			
			with open(outputname, "a") as topsequence:
				topsequence.write(output)

			lencheck = 1


			print j + "  :Buffer too short!"
		else:
			atl = [x[:-1] for x in atlist]
			atlist = atl
			
			# print j
			# print atlist

	
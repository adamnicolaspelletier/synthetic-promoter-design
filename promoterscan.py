#											-*- Mode: Python -*-
#											-*- coding UTF-8 -*-
# scriptname.py
# Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
# Author :  Adam-Nicolas Pelletier
# Last modified On: 

import numpy as np
import itertools
import random
import os
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()
import pandas as pd
from collections import Counter

########################################################################################################################################################
########################################################## USER INPUT ##################################################################################
n = 900   #promoter length



pgl4RE = [NheI, XhoI, HindIII] #list of RE in order on PGL4 vector, for reference purposes. EcoRV was removed: not compatible with Renilla. 



longpwm = 24

startsequence = "AATTATGCTAGCTGCCAACTCCCCCATTAGTGCTCGACTCCACCCGGAGTGACCTTTGACCCCCTACGCCCACGCACTGTCACGTGACCCTGGTGGGGGCAGATAAGCATGCCCAAATAAGGCAAGTAAACAAAAAAAACACAAGGTCACGTGCCCGGAAGCGGAAGTGACAGCTGCTGCGCCGCCATCTTGGAACATTCTGTTCGTATTCCCAAACCCACCAGAGGGATTGCACAATAAACAGATGGTCGCCCAACCACTGGAAGTTTCGGTTTCGTTTGACCTTGAATTTCCGGGAAATGACCTCAAATGACCTTTTCACGTGCTGACAGCTGTCAGGCCGGAAGTGAGTGTGAGTCATCTAACCACAAATACACCTGTCAGCCACACCCAGGCAGCGGGAAGTGGCCTTGGAAATCCCCTAAAATTGCTGAGTCATGAGAGTGCTGATTGGTCCAAAAAGCTTTCTAGGAATTTAATTATTAAGCAAACACACGTGCTAACCACGTGCCCAGGAAGTGCACAGATAAGGAAACCGCAATAAATTGATTGATGGTGCCCACGTGGTGAAAAAGAGGAAGTGAAACTGGGGATTGCATTTAAATGGACCAATCAGCCGTCAGCACAAAAGTAAACAAAGACCCGGAAGTGGATTTCCTGGAATTTCAGCCAATCAGCGCGGCGGGAAATTGGAAAGTTTCGGTTTCGGCACTTCCGGGAAATGATTGCGCAATCCACCCCCCCACCGGAAGTGATGGGGTGATATGCCACGCCCCTTTTTGATTTGCATATTTTGCTGAGTCAGCATAGCAGGAAGTGACTCCGGCCCCGCCCCCTCCCCCTTGGCACGAATTCCTGGAAACCCCCCCCCCGCCCCCGCACAAACCGCAAAGGTAAGTAGGTCACGTGACCCCACCCCCTGACTTTCACTTTCGCCCCGCCCAGTGACGTCACCGTCTAGACACAAACCACAGACTCGAGTGTAAC"
## This 

letters = {0:"A",1:"C",2:"G",3:"T"}  # dictionary of indexes of each nucleotide for matrices
lettersinv = {"A":0,"C":1,"G":2,"T":3}
########################################################################################################################################################
########################################################################################################################################################

def randomseq(length):   
	""" generates a random n length DNA sequence as a numpy array"""
	matrix = np.zeros( (4, length) )
	index = []
	for i in range (length):    
	    index.append([random.randrange(0,4,1), i]) 
	a = np.array(index)  
	matrix[a[:,0], a[:,1]] = 1
	return matrix

def matrixconverter(seqmatrix, letters): 
	"""generates a DNA sequence string based on a NumPy matrix"""
	a = np.transpose(np.nonzero(np.transpose(seqmatrix))).tolist()
	seqstring = ""
	for i in a:
		seqstring += letters[i[1]]
	return seqstring

def matrixmaker(dnastring):   
	""" Generates a numpy array from a DNA sequence string, made from ones and zeros. 2D representation of a DNA sequence. Complements the matrixconverter function""" 
	matrix = np.zeros( (4, len(dnastring)) )
	index = []
	for i in range(len(dnastring)):    
		index.append([lettersinv[list(dnastring)[i]], i]) 
	a = np.array(index)  
	matrix[a[:,0], a[:,1]] = 1
	return matrix



def remove_duplicates(values):
	"""Removes duplicate values from a list"""
	output = []
	seen = set()
	for value in values:
		# If value has not been encountered yet,
		# ... add it to both list and set.
		if value not in seen:
			output.append(value)
			seen.add(value)
	return output

def pwmwalk(pwm, sequence): ### useless for now. could be implemented
	for j in k:	
			alphapos = 0
			alphascore = -1000		
			for i in range(n):
				try:
					betapos = i 
					betascore = np.sum(seqmatrix[:,betapos:(betapos + k[j].shape[1])] * k[j])
					if betascore > alphascore:
						alphascore = betascore
						alphapos = betapos
				except ValueError:
					pass
			if alphascore < topscoredic[j] * 1.25:
				posdic[j] = alphapos 
				scoredic[j] = alphascore - 10000
			else:
				posdic[j] = alphapos
				scoredic[j] = alphascore

def regenerateseq(degeneratestring):
	""" Generates an iterable list of possible sequences in unambiguous IUPAC format based on a IUPACdegenerate sequence. Ex: ATN = [ATT, ATG, ATC, ATA].
	The idea relies on the fact that each IUPAC letters that yields 2 distinct nucleotides has 2 Sequences. So a string with 2 such letters has 2 ^2 sequence. 
	n letters has 2^n sequences. The same logic applies for 3 letters, or 4 letters, with 3^n and 4^n The total number of sequences is thus 2^n(2) * 3^n(3) * 4^n(4). """
	##IUPAC dictionary
	IUPAC1 = ["A", "T", "C", "G"]
	IUPAC2 = ["R", "Y","S","W","K","M"]
	IUPAC3 = ["B", "D", "H", "V"]
	IUPAC4 = "N"

	iupacdna2 = {"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["A","C"]}
	iupacdna3 = {"B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["C","G","T"]}

	twocount = 0
	threecount = 0 
	fourcount = 0
	seq = list(degeneratestring)
	dictio = Counter(seq)
	for i in dictio:
		twocount += IUPAC2.count(i) * dictio[i]
		threecount += IUPAC3.count(i) * dictio[i]
		fourcount += IUPAC4.count(i) * dictio[i]
	possib = (2**twocount) * (3**threecount) * (4**fourcount)  ## calculates the number of possible sequences from the degenerate sequence
	
	seqcomb = []
	for i in seq:
		if IUPAC1.count(i) == 1:
			seqcomb.append(list(i*possib))
			
		elif IUPAC2.count(i) == 1:
			seqcomb.append(list(list((iupacdna2[i][0])*(possib/2))) + list(list((iupacdna2[i][1])*(possib/2))))

		elif IUPAC3.count(i) == 1:
			seqcomb.append(list(list((iupacdna3[i][0])*(possib/3))) + list(list((iupacdna3[i][1])*(possib/3))) + list(list((iupacdna3[i][2])*(possib/3))))
		
		else:
			seqcomb.append(list("A" * (possib/4)) + list("T" * (possib/4)) + list("C" * (possib/4))+ list("G" * (possib/4)))
	seqcombbeta = []
	for i in seqcomb:
		seqcombbeta.append(list(set(i)))

	seqcombdelta = list(itertools.product(*seqcombbeta))
	seqfinal = []
	
	for i in seqcombdelta:
		a = list(i)
		b = matrixmaker("".join(a))
		seqfinal.append(b)
	return seqfinal

def REcheck(seqmatrix,letters,vectorsites):  
	""" Verify if a random generated sequence has pgl4 restriction sites. If it does, return 0. if it doesnt, return 1"""
	seqmat = matrixconverter(seqmatrix, letters)
	pgl4REstr = []
	for i in list(pgl4RE):
		pgl4REstr.append(str(i))

	resites = []
	fseqf = Seq(seqmat, amb)
	resites.append(len(NheI.search(fseqf)))
	resites.append(len(XhoI.search(fseqf)))
	resites.append(len(EcoRV.search(fseqf)))
	resites.append(len(HindIII.search(fseqf)))

	dic = dict(zip(list(pgl4REstr), resites))
	pgl4dict = dict(zip(list(pgl4REstr), list(pgl4RE)))
	checklist = []
	for i in list(itertools.combinations(list(pgl4REstr),2)):
		somme = [dic[list(i)[0]], dic[list(i)[1]]]
		checklist.append(sum(somme))
	if min(checklist) == 0:  # in other words: if there is at least one pair of RE that has zero cleavage sites in the sequence, return 1 to keep that sequence.
		return 1	
	else:
		return 0

def REchoose(seqmatrix,letters,vectorsites):  
	""" Verify if a random generated sequence has pgl4 restriction sites. return all valid RE combinations for a given sequence"""
	seqmat = matrixconverter(seqmatrix, letters)
	pgl4REstr = []
	for i in list(pgl4RE):
		pgl4REstr.append(str(i))

	resites = []
	fseqf = Seq(seqmat, amb)
	resites.append(len(NheI.search(fseqf)))
	resites.append(len(XhoI.search(fseqf)))
	resites.append(len(EcoRV.search(fseqf)))
	resites.append(len(HindIII.search(fseqf)))

	dic = dict(zip(list(pgl4REstr), resites))
	pgl4dict = dict(zip(list(pgl4REstr), list(pgl4RE)))
	checklist = []
	redict = {}
	for i in list(itertools.combinations(list(pgl4REstr),2)):
		somme = [dic[list(i)[0]], dic[list(i)[1]]]
		redict[i] = max(somme)
	
		checklist.append(sum(somme))
	
	return redict	

def numpyreplace(receivepwm,inputpwm, position):
	return np.concatenate((receivepwm[:,:position], inputpwm, receivepwm[:,position + len(inputpwm.transpose()) :] ), axis=1)

def coreseqsort(coremotiflist):
	""" Generates a list of Numpy arrays representing core DNA sequences in a random order. This random order is important to prevent the same TFs being always added first, 
	increasing their likeliness to be overwritten by later TF in the list"""
	random.shuffle(coremotiflist)
	regenlist = []
	coreseqlist = []
	for j in coremotiflist:
		regenlist.append(regenerateseq(j))
	for i in regenlist:
		coreseqlist.append(i[random.randint(0,len(i)-1)])
	return coreseqlist

def promotergenerator(length, coremotiflist,corespec):
	s = randomseq(length)

	for i in coreseqsort(coremotiflist):
		pos = length - len(i.transpose())
		s = numpyreplace(s, i, random.randint(0,pos))
		
	for i in corespec:
		posspec = length + corespec[i]
		regen = regenerateseq(i)
		s = numpyreplace(s, regen[random.randint(0,len(regen)-1)], posspec)
		
	return s
		

### need to integrate regeneratespec into the corespec portion


def main():
	pass



# Generate the list of PWMs to use
filelist = []
for i in os.listdir("Homosapiens_CISBP/pwmsUnif"):
	if i.endswith(".txt"):
		filelist.append("Homosapiens_CISBP/pwmsUnif"+"/"+i)

dflist = []
for i in filelist:
	pwm = np.loadtxt(i, skiprows=1)
	dflist.append(pwm[:,1:].transpose())

k = {}


# Generate a dictionary of filenames as keys, and pwm numpy arrays as values to iterate over
for i in filelist:
	k[i] = np.loadtxt(i, skiprows=1)[:,1:].transpose()



mattop = np.zeros( (4, len(startsequence)) )
index = []
for i in range(len(startsequence)):    
	index.append([lettersinv[list(startsequence)[i]], i]) 
a = np.array(index)  
mattop[a[:,0], a[:,1]] = 1


topposdic = {}
topscoredic = {}


#### Determine the best case scenario score for the most probable promoter as a reference value. Useful to use as an asymptote for the monte carlo iteration . 
for j in k:	
	toppos = 0
	topscore = -1000	
	for i in range(len(startsequence)):
		try:
			candpos = i 
			candscore = np.sum(mattop[:,candpos:(candpos + k[j].shape[1])] * k[j])
			if candscore > topscore:
				topscore = candscore
				toppos = candpos
		except ValueError:
			pass
		topposdic[j] = toppos
		topscoredic[j] = topscore

print sum(topscoredic.values())


pandas = "scorethresholds.txt"
scorethrdf = pd.read_csv(pandas, sep="\t", index_col = "Filename")
scorethrdict = scorethrdf["Score_threshold"].to_dict()


loop = 0


scoretotal = -10000000
seqfinal = ""

omegapos = ""
omegascore = -10000000
omegaposdic = {}
omegascoredic = {}
omegaalign = 0

with open("montecarlolog.txt", "w") as matlablog:
    matlablog.write("0,-600000\n")




### Generate a stable matrix for starting. 
c = matrixmaker(startsequence)



## generate a list of all seeded promoter files inthe output directory to determine which number to add for this output. 
promoterlist = []
for i in os.listdir("promoter_files/seeded/"+str(n)+"bp/"):
	if i.endswith(".txt"):
		promoterlist.append("promoter_files/seeded/"+str(n)+"bp/"+i)

promodict = {} #generate a dictionary of files as keys, and promoter sequences as values. Will be iterated over to find the best promoter. 
for i in promoterlist:
	with open(i,"r") as gg:
		promodict[i] = gg.read()




pversions = [] # lsit of promoter version numbers from the promoterlist
for i in promoterlist:
	pversions.append(int(list(i)[-5]))
f  = max(pversions) + 1 # current version for this promoter run.

del promodict["promoter_files/seeded/900bp/montecarlosequence900bp_1.txt"]
for p in promodict:
	if promodict[p] == "":
		del promodict[p]




finalfile = ""

for p in promodict:
	posdic = {}
	scoredic = {}
	# try:
	seqmatrix = matrixmaker(promodict[p])
	align = 0
	for j in k:	
		alphapos = 0
		alphascore = -1000		
		for i in xrange(n):
			try:
				betapos = i 
				betascore = np.sum(seqmatrix[:,betapos:(betapos + k[j].shape[1])] * k[j])
				if betascore > alphascore:
					alphascore = betascore
					alphapos = betapos
			except ValueError:
				pass
		if alphascore < scorethrdict[j]:
			break
		
		else:
			align += 1
			posdic[j] = alphapos
			scoredic[j] = alphascore
	print p, sum(scoredic.values())
# 		if sum(scoredic.values()) > sum(topscoredic.values()):
# 			pass
# 		#print p + "," + str(sum(scoredic.values()))  #Uncomment if necessary to see distriution of scores for all promoter files.
# 		elif sum(scoredic.values()) > omegascore and align == len(k):
# 			omegascore = sum(scoredic.values())
# 			seqfinal = matrixconverter(seqmatrix, letters)
# 			omegascoredic = scoredic
# 			omegaposdic = posdic
# 			omegaalign = align
# 			finalfile = p

# 		elif align < len(k):
# 			alphaalign = align
# 	except KeyError:
# 		pass
	
	

# 	#output = "%s,%s\n" %(loop, omegascore) 
# 	#with open("montecarlolog.txt", "a") as matlablog:
# 	#	matlablog.write(output)
# print omegascore, finalfile
# print REchoose(matrixmaker(seqfinal), letters, pgl4RE)
		




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
from collections import Counter
import string


########################################################################################################################################################
########################################################## USER INPUT ##################################################################################
## GOAL: Determine the correct oritentation for a list of PWMs based on published CORE sequences. Each PWM and its reverse complement will be assayed, 
## only the PWM with the best score will be kept for the following steps. 


seq = "ATGCTGA"
file1 = "Homosapiens_CISBP/pwmsLOG/M6072_1.01.txt"
file2 = "Homosapiens_CISBP/pwmsLOG/M0374_1.01.txt"

pandas = "Homosapiens_CISBP/pwmsLOG/1.TFDB.txt"

testmotif = "WGATARN"

longpwm = 24
########################################################################################################################################################
########################################################################################################################################################

dic = {0:"A",1:"C",2:"G",3:"T"}  # dictionary of indexes of each nucleotide for matrices
lettersinv = {"A":0,"C":1,"G":2,"T":3}

#pwms = pd.read_csv(pandas, sep="\t", index_col = "DBID")
#pwmmot = pwms.set_index("Motif_ID")
#dictlen = pwmmot["Motif_Length"].to_dict()
#pwmdf = pwms.groupby("TF_Name").apply(lambda tdf: pd.Series(  dict([[vv,tdf[vv].unique().tolist()] for vv in tdf if vv not in ['TF_Name']])  )) 
#pwmdf["maxpwm"] = pd.Series(1, index = pwmdf.index)
#pwmdict = pwmdf['Motif_ID'].to_dict()
#pwmmaxtemp = pwmdf['Motif_Length'].to_dict()


def revcomplement(dnastring):
	a = dnastring
	b = a.replace("A", "X")
	c = b.replace("C", "Y")
	d = c.replace("T", "A")
	e = d.replace("G", "C")
	f = e.replace("X", "T")
	g = f.replace("Y", "G")
	h = g[::-1]
	return h

def randomseq(length):   # generate a random n length DNA sequence as a matrix
	""" Generate a random n length DNA sequence as a numpy array"""
	matrix = np.zeros( (4, length) )
	index = []
	for i in range (length):    
	    index.append([random.randrange(0,4,1), i]) 
	a = np.array(index)  
	matrix[a[:,0], a[:,1]] = 1
	return matrix

def matrixconverter(seqmatrix, dic): 
	""" Generates a DNA squence string based on a numpy array made fro ones and zeros. Complements the matrixmaker function  """
	a = np.transpose(np.nonzero(np.transpose(seqmatrix))).tolist()
	seqstring = ""
	for i in a:
		seqstring += dic[i[1]]
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

def matrevcompl(pwmf):
	"""Takes a pwm in the theoritical Forward orientation (without headers and pos) in one orientation and reverse complements it. Result is numpy array with """
	pwmr = np.fliplr(pwmf)
	pwmrev = np.fliplr(pwmr.transpose())
	return pwmrev.transpose()
	
def topseq(filename):  
	""" Generates a string DNA sequence from a pwm file to represent the most probable DNA sequence representing that pwm. """
	tf1 = np.loadtxt(filename, skiprows=1)
	tf2 = tf1[:,1:].transpose()
	indexes = np.argmax(tf2,axis=0).tolist()  #generates a list of indexes of the highest probability base per position. 
	seqstring = ""
	for i in indexes: #generates a string nucleotide sequence corresponding top nucleotide for each position
		seqstring += dic[i]
	return seqstring # returns the top motif

def topseqrev(filename):  
	""" Generates a string DNA sequence from a pwm file to represent the most probable DNA sequence representing that pwm. """
	tf1 = np.loadtxt(filename, skiprows=1)
	tf2 = matrevcompl(tf1[:,1:].transpose())
	indexes = np.argmax(tf2,axis=0).tolist()  #generates a list of indexes of the highest probability base per position. 
	seqstring = ""
	for i in indexes: #generates a string nucleotide sequence corresponding top nucleotide for each position
		seqstring += dic[i]
	return seqstring # returns the top motif



def remove_duplicates(values):
	""" Takes a list of values and removes duplicates to return a list of values with only a single occurence. Does not delete all occurences of a duplicate, only the extra ones"""
	output = []
	seen = set()
	for value in values:
		# If value has not been encountered yet,
		# ... add it to both list and set.
		if value not in seen:
			output.append(value)
			seen.add(value)
	return output

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

def pwmwalk(pwm, sequence, pos): 
	""" Tests each position of a SHORT DNA sequence as a numpy array against a pwm array. Returns the best alignment score. 
	The Pos argument lets you define over which portion of the pwm you want to align. Add 0 for the whole thing. 
	Useful when you add buffer random sequences on the site to align earlier than the 0 position for both array and sequence, or further than the length of the array"""
	alphapos = 0
	alphascore = -1000		
	for i in xrange(pos,len(pwm.transpose())-pos):
		try:
			betapos = i 
			betascore = np.sum(sequence * pwm[:,betapos:(betapos + sequence.shape[1])])
			if betascore > alphascore:
				alphascore = betascore
				alphapos = betapos
		except ValueError:
			pass
	return [alphascore,alphapos]

def pwmtransfo(array): #this function transforms the pwm in a format (PFM) that can be read and interpreted by the Biopython.mMotif module
	pwm = array
	pwm = pwm * 100
	return pwm

def palindromedetect(sequence):
	pass

# pwm = np.loadtxt(file1, skiprows=1)
# pwmf = pwm[:,1:].transpose()




### Generate a df, dict and list of core motifs from the litterature to align PWMs with. 
coremotifdf = pd.read_csv("TF_core_motifs.txt", sep="\t", index_col = "Gene") 
core = coremotifdf[coremotifdf['Evidence'].isin(["D","I"])]
corelist = core["Motif"].values.tolist()
coredict = core["Motif"].to_dict()
corekeys = coredict.keys()


#generate a list of all files in the original directory. 
filelist = []
for i in os.listdir("Homosapiens_CISBP/pwmsLOGorient"):
	if i.endswith(".txt"):
		filelist.append(i)

filelistclean = []
for i in filelist:
	filelistclean.append(i.replace(".txt",""))


#Import a dataframe of Core motifs from the litterature
core = coremotifdf[coremotifdf['Evidence'].isin(["D","I"])]
#genenames = core["Motif"].to_dict().keys()
coreint = core["Motif"].to_dict()  # dict of Gne names as keys, core motif as values
 



#Import a list of all file vs Gene name associations to filter core motifs against. 
cisbpDB = pd.read_csv('Homosapiens_CISBP/pwmsLOG/1.TFDB.txt', sep="\t", index_col="Motif_ID")  # make a datafram from the list of tf from CisBP
cisbpdict = cisbpDB["TF_Name"].to_dict()   # dictionary with file name as keys, and gene names as values
genes = cisbpdict.values()
genenames = remove_duplicates(genes)
cisbpcore = {}  # Dictionary of file names as keys, and core motifs as values.




coredict = {}
for k in genenames:
	try:
		coredict[k] = coreint[k]
	except KeyError:
		pass


for i in cisbpdict:
	cisbpcore.update({i:coredict[cisbpdict[i]]})



### Iterate all files one by one, and test all possible core motifs for that gene (if they are degenerate) 
mat =  np.ones( (4, 20) ) /4
matlog2 = np.log2((mat/0.25)) 
matlog = matlog2 *2



fcount = 0
rcount = 0
for i in cisbpcore:
	motif = cisbpcore[i]
	motiflist = regenerateseq(motif)
	pwm = np.loadtxt("Homosapiens_CISBP/pwmsLOG/"+i+".txt", skiprows=1)
	pwmf = pwm[:,1:].transpose()
	pwmr = matrevcompl(pwmf)
	iterpwmf = np.concatenate((matlog, pwmf, matlog), axis=1)
	iterpwmr = np.concatenate((matlog, pwmr, matlog), axis=1)
	lenpwm = len(iterpwmf.transpose())


	fscore = -1000
	fpos = 20
	rscore = -1000
	rpos = 20
	for j in regenerateseq(cisbpcore[i]):
		betaf = list(pwmwalk(iterpwmf,j, 17)) 
		betar = list(pwmwalk(iterpwmr,j, 17))

		betafscore = betaf[0]
		betafpos = betaf[1]
		betarscore = betar[0]
		betarpos = betar[1]
	
		if betafscore > fscore and betafpos >= -18 and betafpos <= (lenpwm-18) :
			fscore = betafscore
			fpos = betafpos
		else:
			pass

		if betarscore > rscore and betarpos >= 18 and betarpos <= (lenpwm-18) :
			rscore = betarscore
			rpos= betarpos
		else:
			pass

	filename = 'Homosapiens_CISBP/pwmsLOG/%s.txt' % i
	outfilename = "Homosapiens_CISBP/pwmsLOGorient/%s.txt" % i 

	PFMinfilename = "Homosapiens_CISBP/pwms_all_motifs/%s.txt" % i
	PFMoutfilename = "Homosapiens_CISBP/PFM_transform/%s.pfm" % i

	pfm = np.loadtxt(PFMinfilename, skiprows=1)[:,1:]
	pfm = pfm.transpose()
	pfmf = pwmtransfo(pfm)
	pfmr = matrevcompl(pfmf)

	pwmf = pwm[:,1:].transpose()
	cisbpmat = pd.read_csv(filename, sep="\t", index_col = "Pos") # make a datafram from the list of tf from CisBP
	columns = ["A", "C", "G", "T"]	 
	index = range(1,len(pwmr.transpose())+1)
	cisbpmatr = pd.DataFrame.from_records(matrevcompl(cisbpmat), index = index, columns = columns)
	cisbpmatr.index.names = ["Pos"]

	


	if rscore > fscore:
		cisbpmatr.to_csv(outfilename, sep='\t')
		np.savetxt(PFMoutfilename,pfmr, delimiter=" " )
		print "F: "+ topseq(filename)
		print "R: "+ topseqrev(filename) + "   Y" 
		print "Motif: " + motif + "\n" + "\n"

	else: 
		cisbpmat.to_csv(outfilename, sep='\t')
		np.savetxt(PFMoutfilename,pfmf, delimiter=" " )
		print "F: "+ topseq(filename) + "   Y" 
		print "R: "+ topseqrev(filename)
		print "Motif: " + motif + "\n" + "\n"
	#print str(i) + "\n" + str(fscore) + "\t" + str(fpos) + "\t" + str(lenpwm-40) + "\t" + str(len(motif)) + "\n" + str(rscore) + "\t" + str(rpos) + "\t" + str(lenpwm-40) + "\t" + str(len(motif))


#print cisbpdict
#print cisbpclean
#dflist = []
#for i in filelist:
#	pwm = np.loadtxt(i, skiprows=1)
#	dflist.append(pwm[:,1:].transpose())




##
##Idea: Generate a list of all files, lookup the corresponding TF. Test the score of each substring core motif against that pwm and its reverse complement
#, keep the highest score (if it aligns between 20 and 20+ len(pwm.tranpose()). The highest score between the reverse and correct orientation will be kept. if the scores are equal, keep the forward orientation (unlikely)

# 
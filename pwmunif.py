#											-*- Mode: Python -*-
#											-*- coding UTF-8 -*-
# pwmunif.py
# Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
# Author :  Adam-Nicolas Pelletier
# Last modified On: 2015-02-16

import numpy as np
import pandas as pd
import itertools
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import math
from Bio import motifs
from collections import Counter

########################################################################################################################################################
########################################################## USER INPUT ##################################################################################

pandas = "Homosapiens_CISBP/pwmsLOG/1.TFDB.txt"
gataca = ["M2192_1.01","M4351_1.01","M4385_1.01","M4420_1.01","M4449_1.01","M4508_1.01","M6350_1.01"]

########################################################################################################################################################
########################################################################################################################################################

lettersinv = {"A":0,"C":1,"G":2,"T":3}




def montecarlomatrix(filenames):
	num = range(1, len(filenames) + 1)
	return num

def matrixmaker(dnastring):   
	""" Generates a numpy array from a DNA sequence string, made from ones and zeros. 2D representation of a DNA sequence. Complements the matrixconverter function""" 
	matrix = np.zeros( (4, len(dnastring)) )
	index = []
	for i in range(len(dnastring)):    
		index.append([lettersinv[list(dnastring)[i]], i]) 
	a = np.array(index)  
	matrix[a[:,0], a[:,1]] = 1
	return matrix

def mastermatrix(arrayname):
	backgrndlog = math.log(0.25/0.25)
	matrix = np.ones( (50, 4) ) * (backgrndlog) # where the latter is the log-likelihood of having a random base (0.25) at a given position
	pwm = (np.loadtxt(arrayname, skiprows=1))[:,1:]
	matrix[20:len(pwm)+ 20] += (pwm - backgrndlog)
	return matrix

def itermatrix(arrayname, position):
	backgrndlog = math.log(0.90/0.25)
	matrix = np.ones( (50, 4) ) * backgrndlog
	pwm = (np.loadtxt(arrayname, skiprows=1))[:,1:]
	matrix[position:len(pwm)+ position] += (pwm - backgrndlog)
	return matrix

def scorealign(array, pwmlist):
	matrixlist = pwmlist
	#matrixlist.remove(array)
	c = mastermatrix(array)
	poslist = []

	scorealpha = 0
	posalpha = 0
	scoreomega = 0
	posomega = 0
	
	for i in matrixlist:
		for j in range(30):
			d = itermatrix(i, j)
			scoreomega = np.sum(c * d)
			posomega = j
			if int(scorealpha) < int(scoreomega):
				scorealpha = scoreomega
				posalpha = posomega
		poslist.append([array,i,posalpha - 20, len(array)])
	return poslist


def pwmtransfo(matrix): #this function transforms the pwm in a format that can be read and interpreted by the Biopython.mMotif module
	pwm = (np.loadtxt(matrix, skiprows=1))[:,1:]
	pwm = pwm.transpose()
	pwm = pwm * 100
	return pwm

def pearsonpwm(pwm1,pwm2): #this function computes the Pearson coefficient between 2 pwm
	cisbpmat = motifs.read(open(pwm1), "pfm")
	cisbpmat.pseudocounts = 3.0
	pwmnp1 = (np.loadtxt(pwm1, skiprows=1))
	tf1 = cisbpmat.pssm
	cisbpmat2 = motifs.read(open(pwm2), "pfm")
	cisbpmat2.pseudocounts = 3.0
	pwmnp2 = (np.loadtxt(pwm2, skiprows=1))
	tf2 = cisbpmat2.pssm
	distance, offset = tf2.dist_pearson(tf1)
	return [pwm1, 1 - distance , math.fabs(offset) , len(np.transpose(pwmnp1))], [pwm2, 1 - distance , math.fabs(offset), len(np.transpose(pwmnp2))]




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

def corescan(filename, core):
	""" Takes a core consensus motif and scans a numpy PWM for the best alignment site, then returns the score. """
	mat =  np.ones( (4, 20) ) /4 #Creates a buffer array to flank the main PWM for iterationbs over longer distances. 
	matlog2 = np.log2((mat/0.25)) 
	matlog = matlog2 *2  # This ensures that this buffer sequence is less likely to be aligned with the sequence to prevent unspecific interactions, yet still allows alignment.

	pwm = np.loadtxt(filename, skiprows=1)
	pwmf = pwm[:,1:].transpose()

	iterpwm = np.concatenate((matlog, pwmf, matlog), axis=1) #iterable PWM , flanked by buffer arrays
	lenpwm = len(iterpwm.transpose())

	score = -1000
	pos = 20
	for j in regenerateseq(core):
		beta = list(pwmwalk(iterpwm,j, 17)) 

		betascore = beta[0]
		betapos = beta[1]
		
		if betascore > score and betapos >= -18 and betapos <= (lenpwm-18) :
			score = betascore
			pos = betapos
		else:
			pass

	return score

def topwm(filenames, pwmlist, geneID, core):   
	""" this function retrieves the filenames of the pwms with A. The highest PEarson coefficients, and then compares them to find the one with the second 
	highest coefficient. This allows to select the top pwm. """
	names= ["PWM","Pearson_coeff","Offset", "Motif_Length"]
	filelist = []  #input filelist
	ofilelist = []  #output filelist
	inputlist = filenames
	for i in pwmlist:
		filelist.append("Homosapiens_CISBP/pwms_all_motifs/%s.txt" % i) 
		ofilelist.append("Homosapiens_CISBP/PFM_transform/%s.pfm" % i) 

	d = dict(zip(filelist,ofilelist))
	#dend = dict(zip(ofilelist,filelist))
	dmotif = dict(zip(ofilelist,filenames))
	
	pears = []
	for i in list(itertools.combinations(ofilelist,2)):
		pears.append(pearsonpwm(i[0],i[1])[0])
		pears.append(pearsonpwm(i[0],i[1])[1])

	df = pd.DataFrame(pears)
	df.columns = names
	df = df[(df["Offset"] < 10)]
	df1 = df[(df["Pearson_coeff"] == df["Pearson_coeff"].max())]
	df2 = df.groupby("PWM").apply(lambda tdf: pd.Series(  dict([[vv,tdf[vv].unique().tolist()] for vv in tdf if vv not in ["PWM"]])  )) 

	dictpears = df2['Pearson_coeff'].to_dict()
	dflist = df1["PWM"].tolist()
	
	for key in dictpears:
		dictpears[key].sort()

	filtz = {k:v for (k,v) in dictpears.iteritems() if dflist[0] in k or dflist[1] in k}
	try: 

		pwm1 = pd.read_csv(dmotif[dflist[0]], sep="\t", index_col = "Pos")
		pwm2 = pd.read_csv(dmotif[dflist[1]], sep="\t", index_col = "Pos")
		if len(pwm1) < len(pwm2):  #takes the second last value to determine which pwm in the pair ranks best, since the top value is the same for both. 
			return dmotif[dflist[1]]
		else:
			return dmotif[dflist[0]]  # uncomment the last section to score based on core motif similarity instead. 



		# if corescan(dmotif[dflist[0]], core) < corescan(dmotif[dflist[1]], core):  #takes the second last value to determine which pwm in the pair ranks best, since the top value is the same for both. 
		# 	return dmotif[dflist[1]]
		# else:
		# 	return dmotif[dflist[0]]
	except IndexError:
		return dmotif[dflist[0]]

def unifarray(pwmlist, geneID, dictlen, core):   
	"""creates a unified array busing all pwm from a common TF. The function imports all the numpy arrays from a list of files for a specified geneID, 
	aligns them with one another to extract the Pearson coefficients for all PWM pairs. It then keeps the highest pair. 
	To decide which PWM from that pair is kept, it aligns each PWM with the CORE consensus sequence to ultimately have the PWM that corresponds best both to all other
	 similar experiemnts and to the consensus litterature"""
	gene = geneID
	inputlist= []
	motiflist = pwmlist
	for i in pwmlist:
		app = "Homosapiens_CISBP/pwmsLOGorient/%s.txt" % i
		inputlist.append(app) 
	outfilename = "Homosapiens_CISBP/pwmsUnif/%s.txt" % gene
	if len(inputlist) == 1: #if there is only one PWM for a given gene: make that PWM the final pwm and give it the genes name
		pwm = pd.read_csv(inputlist[0], sep="\t", index_col = "Pos") 
		output = "%s,%s\n" % (gene,inputlist[0])
		with open("pwmuniflog.txt", "a") as pwmuniflog:
				pwmuniflog.write(output) 
		return pwm.to_csv(outfilename, sep='\t')

	elif len(inputlist) == 2: #if there are 2 pwm: pearson coefficient cannot help, since it will be the same. take the one that aligns best with the CORE sequence. 
		# if corescan(inputlist[0], core) > corescan(inputlist[1], core):
		# 	pwm = pd.read_csv(inputlist[0], sep="\t", index_col = "Pos") 
		# 	output = "%s,%s\n" % (gene,inputlist[0])
		# 	with open("pwmuniflog.txt", "a") as pwmuniflog:
		# 		pwmuniflog.write(output)
		# 	return pwm.to_csv(outfilename, sep='\t')
		# elif corescan(inputlist[0], core) < corescan(inputlist[1], core):
		# 	pwm = pd.read_csv(inputlist[1], sep="\t", index_col = "Pos") 
		# 	output = "%s,%s\n" % (gene,inputlist[1])
		# 	with open("pwmuniflog.txt", "a") as pwmuniflog:
		# 		pwmuniflog.write(output)
		# 	return pwm.to_csv(outfilename, sep='\t')
		# else:  #If both align equally, take the first one and move on. 
		# 	pwm = pd.read_csv(inputlist[0], sep="\t", index_col = "Pos") 
		# 	output = "%s,%s\n" % (gene,inputlist[0])
		# 	with open("pwmuniflog.txt", "a") as pwmuniflog:
		# 		pwmuniflog.write(output)
		# 	return pwm.to_csv(outfilename, sep='\t')
		#print len(pd.read_csv(inputlist[0], sep="\t", index_col = "Pos"))


####################### Comment after this line to remove size choice and uncomment previous block to use core motif selection

		pwmy = pd.read_csv(inputlist[0], sep="\t", index_col = "Pos")
		pwmz = pd.read_csv(inputlist[1], sep="\t", index_col = "Pos")
		if len(pwmy) >= len(pwmz):
			pwm = pd.read_csv(inputlist[0], sep="\t", index_col = "Pos") 
			output = "%s,%s\n" % (gene,inputlist[0])
			with open("pwmuniflog.txt", "a") as pwmuniflog:
				pwmuniflog.write(output)
			return pwm.to_csv(outfilename, sep='\t')
		else:
			pwm = pd.read_csv(inputlist[1], sep="\t", index_col = "Pos") 
			output = "%s,%s\n" % (gene,inputlist[1])
			with open("pwmuniflog.txt", "a") as pwmuniflog:
				pwmuniflog.write(output)
			return pwm.to_csv(outfilename, sep='\t')
		
	else:
		motif = topwm(inputlist, pwmlist, geneID, core)  
		name = motif
		pwm = pd.read_csv(name, sep="\t", index_col = "Pos") 
		output = "%s,%s\n" % (gene,name)
		with open("pwmuniflog.txt", "a") as pwmuniflog:
			pwmuniflog.write(output)
		return pwm.to_csv(outfilename, sep='\t')


motiflist = []

 
pwms = pd.read_csv(pandas, sep="\t", index_col = "DBID")
pwmmot = pwms.set_index("Motif_ID")
pwms = pwms[pwms["Motif_Length"] <= 25]  # this is necessary to prevent using VERY long pwms. the value can be reduced as well. Otherwise... memory error on scoredistrib. 
dictlen = pwmmot["Motif_Length"].to_dict()



 
### Generate a df, dict and list of core motifs from the litterature to align PWMs with. 
coremotifdf = pd.read_csv("TF_core_motifs.txt", sep="\t", index_col = "Gene")   
core = coremotifdf[coremotifdf['Evidence'].isin(["D","I"])]
coredict = core["Motif"].to_dict()


pwmsblank = pd.DataFrame([["BLANK", "BLANK", 0, 0]], columns=["DBID", "TF_Name", "Motif_ID", "Motif_Length"], )
pwmsblank = pwmsblank.set_index("DBID")



for i in coredict:
	try:
		pwmstemp = pwms[pwms["TF_Name"] == i]
		if len(pwmstemp) > 1:
			pwmstemp2 = pwmstemp[pwmstemp["Motif_Length"] < (len(coredict[i]) * 2.5)]
			pwmsfilt = pwmsblank.append(pwmstemp2)
			pwmsblank = pwmsfilt
			if len(pwmstemp2) < 1:
				pwmstemp2 = pwmstemp[pwmstemp["Motif_Length"] <= 20]
				pwmsfilt = pwmsblank.append(pwmstemp2)
				pwmsblank = pwmsfilt
			
			#print i, pwmstemp2
		elif len(pwmstemp) == 1: 
			pwmsfilt = pwmsblank.append(pwmstemp)
			pwmsblank = pwmsfilt
		

	except KeyError:
		pass

pwmsfilt = pwmsblank[pwmsblank["TF_Name"] != "BLANK"]


pwmdf = pwmsfilt.groupby("TF_Name").apply(lambda tdf: pd.Series(  dict([[vv,tdf[vv].unique().tolist()] for vv in tdf if vv not in ['TF_Name']])  )) 
pwmdf["maxpwm"] = pd.Series(1, index = pwmdf.index)
pwmdict = pwmdf['Motif_ID'].to_dict()
pwmmaxtemp = pwmdf['Motif_Length'].to_dict()

print pwmdf



keycount = 0
with open("pwmuniflog.txt", "w") as pwmuniflog:
    pwmuniflog.write("0,0\n")

for key in pwmdict:
	keycount += 1
	print keycount, key
	unifarray(pwmdict[key],key,dictlen,coredict[key])

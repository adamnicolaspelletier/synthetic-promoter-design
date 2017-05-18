import numpy as np
from os import listdir
from os.path import isfile, join
import os.path
import pandas as pd
import itertools

# print listdir("Homosapiens_CISBP/pwms_all_motifs")




###############################################################################################################################
###############################################################################################################################
## User Input


fasta = open("tfcds_oneline_restsites.txt")  #name of the input fasta one liner file. 
fastalist = fasta.readlines()
del fastalist[0] #  "#" This next line if your fasta file has no empty first line. Fastaconvertscript makes an empty line by default. 

###############################################################################################################################
###############################################################################################################################

cisbpDB = pd.read_csv('Homosapiens_CISBP/TF_Information_all_motifs.txt', sep="\t", index_col="DBID")  # make a datafram from the list of tf from CisBP
cisbpclean = cisbpDB[["TF_Name","Motif_ID"]]



a = ""
x = []

for i in fastalist:
	if ">" in i:
		a += i
		a = a.replace("\n","")
		x.append(a) #creates a list of ID
		a = ""


name = ""
ID = []
x.pop()
for i in x:
	name = str(i)
	idsplit = name.split("|")
	ID.append(idsplit[2])
ID.append("SP1")



cisbpclean = cisbpclean[(cisbpclean.Motif_ID != ".") & (cisbpclean["TF_Name"].isin(ID))]


motiflist = []
for index, row in cisbpclean.iterrows():
	motiflist.append(row['Motif_ID'])



motiflen = [] #lsit for length of each motif

for i in motiflist:
	filename = 'Homosapiens_CISBP/pwms_all_motifs/%s.txt' % i
	if os.path.isfile(filename) == True:
		cisbpmat = pd.read_csv(filename, sep="\t", index_col = "Pos") # make a datafram from the list of tf from CisBP
		motiflen.append(len(cisbpmat))
	else:
		motiflen.append(".")


cisbpclean["Motif_Length"] = motiflen


cisbpclean = cis
bpclean[(cisbpclean.Motif_Length > 0)]
print cisbpclean
print cisbpclean[["Motif_Length"]].mean(axis=0)

motiflist = []
for index, row in cisbpclean.iterrows():
	motiflist.append(row['Motif_ID'])



def pwm_norm(dataframe): ## Bayesian uniform prior normalization of PWM to remove zeros. 
	alpha = 0.99
	pp = (1 - alpha) / (4 * (1-0.25)) 
	dataframepr = (dataframe * alpha) + ((1 - 0.25) * pp)
	return dataframepr

def seqloglike(seq, matrix): # log-likelihood calculation for every possible nucleotide sequence of length n
	pos = 1
	loglike = 0
	for i in seq:
		letter = i
		loglike += matrix.loc[pos,letter]
		pos += 1
	return loglike

def pwmtransfo(matrix): # no longer used in script
	"""this function transforms the pwm in a format that can be read and interpreted by the Biopython.mMotif module"""
	pwm = (np.loadtxt(matrix, skiprows=1))[:,1:]
	pwm = pwm.transpose()
	pwm = pwm * 100
	return pwm

#seqoutput = pd.DataFrame(np.random.randn(6,4),index=dates,columns=list('ABCD'))

for i in motiflist:
	motifname = i
	filename = 'Homosapiens_CISBP/pwms_all_motifs/%s.txt' % i
	outfilename = "Homosapiens_CISBP/pwmsselect/%s.txt" % i 
	cisbpmat = pd.read_csv(filename, sep="\t", index_col = "Pos") # make a datafram from the list of tf from CisBP
	cisbpmat = pwm_norm(cisbpmat) 
	cisbpmat.to_csv(outfilename, sep='\t')

for i in motiflist:
	motifname = i
	filename = 'Homosapiens_CISBP/pwms_all_motifs/%s.txt' % i
	outfilename = "Homosapiens_CISBP/pwmsLOG/%s.txt" % i 
	cisbpmat = pd.read_csv(filename, sep="\t", index_col = "Pos") # make a datafram from the list of tf from CisBP
	cisbpmat = pwm_norm(cisbpmat)
	logcisbp = np.log2((cisbpmat/0.25))  
	logcisbp.to_csv(outfilename, sep='\t')

cisbpclean.to_csv("Homosapiens_CISBP/pwmsLOG/1.TFDB.txt", sep='\t')


#	a = ["".join(i) for i in itertools.product("ATCG", repeat=len(logcisbp))]

#	for m in a:
#		seq = m
#		seqlist = list(seq)
#		value = seqloglike(seqlist,cisbpmat)
#		df = 




#IDEA = if value is above threshold x (0.6) ; keeep that letter in  Uppercase. If no value is above 0.6, keep the highest in lower case
#IDEA 2 = if there is no top contender, generate 2 distinct consensus motifs for the top 2 contenders, both in lower case. THAT lower case letter could potentially be used to mix motifs among factors later on.  


import pandas as pd
from Bio.Restriction import *
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()

###############################################################################################################################
###############################################################################################################################
## User Input


fasta = open("tfcds_oneline_test.txt")  #name of the input fasta one liner file. 
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

print cisbpclean["TF_Name"]




#1. Open Cisbp list of TF as a dataframe
#	Slice the Ensembl ID, Gene symbol and motif ID columns. 
#	Make a string out of those 3, tab-delimited. 
#	Make a list of those strings for all rows of the dataframe.
#Loop through that list : 
#	for each item in list:
#		split by tab in a list of 3.
#		open the matrix corresponding to the motif ID
#		retrieve sequence from matrix in IUPAC nomenclature. 
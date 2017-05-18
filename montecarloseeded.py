

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### scriptname.py
### Copyright 2015 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 

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

seq = "ATGCTGA"
file1 = "Homosapiens_CISBP/pwmsLOG/M0186_1.01.txt"
file2 = "Homosapiens_CISBP/pwmsLOG/M0374_1.01.txt"

pgl4RE = [NheI, XhoI, HindIII] #list of RE in order on PGL4 vector, for reference purposes. EcoRV was removed: not compatible with Renilla. 


fivepr = "AATTAT"  # write the sequence you want to add in the 5' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
threepr = "TGTAAC" # write the sequence you want to add in the 3' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended

longpwm = 24

startsequence = "TGGTGACTCATCCTCTTCCCGGAAACATTTCCCGGAAATGAAAAAAACACAATAACCACAAAGGCGGGGTCACGTGACATTTGCATATTGGACCAATCAGCACTCTCTACGCCCACGCACTTGACAGCTGTCATGTTTGCTTAAGCCACGCCCAGAACAGATGGTCACCGGAAGTGAAATTCCAGGAAATCTGTTTATTTACTTACCGAAAAGTGAAACCGGTCACGTGCCCGAAACCGAAACTTAGCAGGAAGTGACTGGGAACATTATGTACCCCAAAAAGGGGCGTGGCATTCACGTGCCATTAATCAACACCAGGGGGCGCCTAGCAACAGATGGCAGCAGCTGTTCATTCATTCATAAATGGACCAATCAGCATGCCCAAATAAGGCAAGCCCCGCCCATTTCCAGGAATTGCAGATAACACTGACGTCACGGCGGCCGCCATCTTGACAGGTGTAGTAAACATGTTTACTTTAATTAAGCGGAAGTGACTTGGGAATACCGGCCGGAAGTGACCCGGAAGTGAGATAACCACTTGCTGACTCAGCAAATTCAAATTTCCCGCCAAATGCAATCCCCCCCCCCCCCCGCCCCCGCACCCCAGGAAGTGCAGCCACACCCAGGCAAAAAAGCGGAAGTAGCGGGAAGTGGCAGCCAATCAGCGCAGGGGATTTCCAAGGACAAACCACAGGGGGGGGTGGCTCCGGGACAGATAAGAACTTTGTTTACTTTTGACGAAACCGAAACTTAAATTGATTGATGGCCTCACCCCACCCCAAACCACAGAAAAATTGCTGAGTCATGGAGGTCATTCAGCACCACGGACAGCGCCAAAAGCTTTCTAGGAAGGGGTCAAAGGTCATGCTGACGGGGTCACGTGACCCTTGTCAATTGCGCAATTTCAAGGTCATTATTGTGCAATTAAGCAGTTCCTCTTGGCACGTGCCCACGTGGTGCGTCTAGACAAACCACAAACCCCA"
## This 

letters = {0:"A",1:"C",2:"G",3:"T"}  # dictionary of indexes of each nucleotide for matrices
lettersinv = {"A":0,"C":1,"G":2,"T":3}
########################################################################################################################################################
########################################################################################################################################################

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

def randomseq(length, format):   
    """ generates a random n length DNA sequence as a numpy array"""
    matrix = np.zeros( (4, length) )
    index = []
    for i in range (length):    
        index.append([random.randrange(0,4,1), i]) 
    a = np.array(index)  
    matrix[a[:,0], a[:,1]] = 1
    if format == "numpy":
        return matrix
    elif format == "string":
        return matrixmaker(matrix)



def weighted_choice_sub(weights):  
    """Returns the index of an element in a list. Allows to use random with weights. Eli bendersky code Copyright.""" 
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i


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

def pwmwalk(pwm, sequence, pos): 
    """ Tests each position of a SHORT DNA sequence as a numpy array against a pwm array. Returns the best alignment score. 
    The Pos argument lets you define over which portion of the pwm you want to align. Add 0 for the whole thing. 
    Useful when you add buffer random sequences on the site to align earlier than the 0 position for both array and sequence, or further than the length of the array"""
    
    alphapos = 0
    alphascore = -1000              
    for i in xrange(pos,len(pwm.transpose())-pos):
        try:
            betapos = i 
            betascore = np.sum(sequence * pwm[:,betapos:(betapos + len(sequence.transpose()))])
            if betascore > alphascore:
                alphascore = betascore
                alphapos = betapos
        except ValueError:
            pass
    return [alphascore,pos]



def regenerateseq(degeneratestring, format):
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
            seqcomb.append(list(i))
            
        elif IUPAC2.count(i) == 1:
            seqcomb.append(list(list((iupacdna2[i][0]))) + list(list((iupacdna2[i][1]))))

        elif IUPAC3.count(i) == 1:
            seqcomb.append(list(list((iupacdna3[i][0]))) + list(list((iupacdna3[i][1]))) + list(list((iupacdna3[i][2]))))
        
        else:
            seqcomb.append(list("A") + list("T") + list("C")+ list("G"))
    seqcombbeta = []
    for i in seqcomb:
        seqcombbeta.append(list(set(i)))

    seqcombdelta = list(itertools.product(*seqcombbeta))
    
    seqfinal = []
    
    if format == "numpy":
        for i in seqcombdelta:
            b = matrixmaker("".join(i))
            seqfinal.append(b)
        return seqfinal

    
    elif format == "string":
        for i in seqcombdelta:
            b = "".join(i)
            seqfinal.append(b)
        return seqfinal

def topseq(numpyarray):  
    """ Generates the top motif from the numpy array based on the position and length of the entry core motif"""

    tf2 = numpyarray
    indexes = np.argmax(tf2,axis=0).tolist()  #generates a list of indexes of the highest probability base per position. 
    seqstring = ""
    for i in indexes: #generates a string nucleotide sequence corresponding top nucleotide for each position
        seqstring += letters[i]

    return seqstring#[alphapos:alphapos+alphalen:] # returns the top motif

def worstseq(numpyarray):  
    """ Generates the top motif from the numpy array based on the position and length of the entry core motif"""

    tf2 = numpyarray
    indexes = np.argmin(tf2,axis=0).tolist()  #generates a list of indexes of the lowest probability base per position. 
    seqstring = ""
    for i in indexes: #generates a string nucleotide sequence corresponding worse   nucleotide for each position
        seqstring += letters[i]

    return seqstring#[alphapos:alphapos+alphalen:] # returns the top motif

def secondseq(numpyarray):  
    """ Generates the 2nd top motif from the numpy array based on the position and length of the entry core motif"""

    tf2 = numpyarray


    #2. 2nd best nucleotide at every position
    tfmax2 = tf2.copy()
    matrix2 = np.zeros( (4, len(tf2.transpose())) ) 
    indmax2 = []
    for i in xrange(len(tf2.transpose())):
        indmax2.append([np.argmax(tfmax2,axis=0).tolist()[i], i])
    mask2 = np.array(indmax2)
    
    tfmax2[mask2[:,0], mask2[:,1]] = -20
    
    index2 = np.argmax(tfmax2,axis=0).tolist() 
    seqstring2 = ""
    for i in index2:
        seqstring2 += letters[i]
    
    return seqstring2



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


def seqreplace(receiveseq,inputseq, position, inputformat):
    """Replaces any portion of a DNA sequence, either in numpy array or string format, into another, at the desired position"""
    if inputformat == "numpy":
        return np.concatenate((receiveseq[:,:position], inputseq, receiveseq[:,position + len(inputseq.transpose()) :] ), axis=1)
    elif inputformat == "string":
        return receiveseq[0:position:] + inputseq + receiveseq[position+len(inputseq)::]


def coreseqsort (coremotiflist, format):
    """ Generates a list of string  core DNA sequences in a random order. This random order is important to prevent the same TFs being always added first, 
    increasing their likeliness to be overwritten by later TF in the list"""
    random.shuffle(coremotiflist)
    coremotiflist.sort(key=len)
    coremotiflist = coremotiflist + endlist
    
   
    regenlist = []
    coreseqlist = []

    if format == "numpy":
        for i in coremotiflist:
            coreseqlist.append(matrixmaker(i)) ### this is wrong
        return coreseqlist

    elif format == "string":
        for i in coremotiflist:
            coreseqlist.append(i)
        return coreseqlist

def promotergenerator(length, coremotiflist,corespec, format):
    s = matrixconverter(randomseq(length, "numpy"),letters)
    for i in coreseqsort(coremotiflist, "string"):
        pos = random.randint(0,length - len(i))
        s = seqreplace(s, i, pos, "string")

    for i in corespec:
        posspec = length + corespec[i]
        regen = regenerateseq(i, "string")
        s = seqreplace(s, regen[random.randint(0,len(regen)-1)], posspec,  "string")
    if format == "string":
        return s
    elif format == "numpy":
        return matrixmaker(s)

def corescan(gene, core):
    """ Takes a core consensus motif and scans a numpy PWM for the best alignment site, then returns the topmotif. """
    mat =  np.ones( (4, 20) ) /4 #Creates a buffer array to flank the main PWM for iterationbs over longer distances. 
    matlog2 = np.log2((mat/0.25)) 
    matlog = matlog2 *2  # This ensures that this buffer sequence is less likely to be aligned with the sequence to prevent unspecific interactions, yet still allows alignment.

    pwmf = kgene[gene]
    

    iterpwm = np.concatenate((matlog, pwmf, matlog), axis=1) #iterable PWM , flanked by buffer arrays
    lenpwm = len(iterpwm.transpose())

    score = -1000
    pos = 20
    topmotif = ""
    for j in regenerateseq(core, "numpy"):
        
        beta = list(pwmwalk(iterpwm,j, 17)) 

        betascore = beta[0]
        betapos = beta[1]
        
        if betascore > score and betapos >= -18 and betapos <= (lenpwm-18) :
            score = betascore
            pos = betapos
            topmotif = j
        else:
            pass

    return topmotif


def corescanupdate(gene, core):
    """ Takes a core consensus motif and scans a numpy PWM for the best alignment site, then returns the topmotif. """
    mat =  np.ones( (4, 20) ) /4 #Creates a buffer array to flank the main PWM for iterationbs over longer distances. 
    matlog2 = np.log2((mat/0.25)) 
    matlog = matlog2 *2  # This ensures that this buffer sequence is less likely to be aligned with the sequence to prevent unspecific interactions, yet still allows alignment.

    pwmf = kgene[gene]
    

    iterpwm = np.concatenate((matlog, pwmf, matlog), axis=1) #iterable PWM , flanked by buffer arrays
    lenpwm = len(iterpwm.transpose())

    score = -1000
    pos = 20
    topmotif = ""
    for j in regenerateseq(core, "numpy"):
        
        beta = list(pwmwalk(iterpwm,j, 17)) 

        betascore = beta[0]
        betapos = beta[1]
        
        if betascore > score and betapos >= -18 and betapos <= (lenpwm-18) :
            score = betascore
            pos = betapos
            topmotif = j
        else:
            pass

    if pos >= 20:
        print topseq(iterpwm[pos:pos+len(topmotif)],matrixmaker(topmotif))


### need to integrate regeneratespec into the corespec portion


def montecarlo():
    pass


#weighted_choices = [2,3,3,2]   # Add weights here. Disabled for now, analysis of our current PWM suggest all bases are represented equally
bases = ["A", "C", "G", "T"]


### Generate a df, dict and list of core motifs from the litterature to align PWMs with. 
#coremotifdf = pd.read_csv("TF_core_motifs.txt", sep="\t", index_col = "Gene") 
coremotifdf = pd.read_csv("TF_core_update.txt", sep="\t", index_col = "Gene") 
core = coremotifdf[coremotifdf['Evidence'].isin(["D","I"])]

coredict = core["Motif"].to_dict()  # dict of Gne names as keys, core motif as values
corelistint = coredict.values()



corefinal = {}
for i in coredict:
    corefinal[i] = regenerateseq(coredict[i], "string")

corespec = {}
corespec[coredict["SP1"]] = -136
corespec[coredict["MAX"]] = -6
corespec[coredict["ZFP36L1"]] = -37
corespec["TATAAA"] = -26

del coredict["ZFP36L1"] #This TF cannot be in the list of random core sequences to be positioned, as it has to be at a specific position. Will be added later.
del coredict["SP1"] #This TF cannot be in the list of random core sequences to be positioned, as it has to be at a specific position. Will be added later.
del coredict["MAX"] #This TF cannot be in the list of random core sequences to be positioned, as it has to be at a specific position. Will be added later.


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
kgene = {}

# Generate a dictionary of filenames as keys, and pwm numpy arrays as values to iterate over. and a dictionary of gene names as keys, with numpy arrays to iterate over. 
for i in filelist:
    gene = str.replace(i,"Homosapiens_CISBP/pwmsUnif/","")    
    gene = str.replace(gene,".txt","")
    k[i] = np.loadtxt(i, skiprows=1)[:,1:].transpose()
    kgene[gene] = np.loadtxt(i, skiprows=1)[:,1:].transpose()



genes = kgene.keys()

# coredictint = {}
# for i in genes:
#     try:
#         print kgene[i]
#         coredictint[i] = corescanupdate(i, coredict[i])
#     except KeyError:
#         pass
        
# print coredictint


corelist = []
endlist = []
for i in coredict:
    try:
        if i == "HOXA6":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        elif i == "DDIT3":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        elif i == "NFE2":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        elif i == "CEBPG":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        elif i == "CEBPA":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        
        elif i == "POU2F2":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        elif i == "KLF7":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        elif i == "MAFF":
            pou2f2 = matrixconverter(corescan(i,coredict[i]),letters)
            endlist.append(pou2f2)
        else:    
            corelist.append(matrixconverter(corescan(i,coredict[i]),letters))
    except KeyError:
        pass


endlist.sort(key = lambda s: len(s))










mattop = np.zeros( (4, len(startsequence)) )
index = []
for i in range(len(startsequence)):    
    index.append([lettersinv[list(startsequence)[i]], i]) 
a = np.array(index)  
mattop[a[:,0], a[:,1]] = 1


topposdic = {}
topscoredic = {}
worstscoredic = {}
secondscoredic = {}

faildic = {}

testno = 0
genelist = []


topsequence = ""
worstsequence = ""
secondsequence = ""
for j in k: 
    topsequence += topseq(k[j])
    worstsequence += worstseq(k[j])
    secondsequence += secondseq(k[j])

topseqmatrix = matrixmaker(topsequence)
worstseqmatrix = matrixmaker(worstsequence)
secondseqmatrix = matrixmaker(secondsequence)


for j in k:    
    alphapos = 0
    alphascore = -1000    


    for i in xrange(len(topsequence)):
        try:
            betapos = i 
            betascore = np.sum(topseqmatrix[:,betapos:(betapos + k[j].shape[1])] * k[j])
            if betascore > alphascore:
                alphascore = betascore
                alphapos = betapos
    
        except ValueError:
            pass
    topscoredic[j] = alphascore


    for i in xrange(len(worstsequence)):
        try:
            betapos = i 
            betascore = np.sum(worstseqmatrix[:,betapos:(betapos + k[j].shape[1])] * k[j])
            if betascore < alphascore:
                alphascore = betascore
                alphapos = betapos
    
        except ValueError:
            pass
    worstscoredic[j] = alphascore

    for i in xrange(len(secondsequence)):
        try:
            betapos = i 
            betascore = np.sum(secondseqmatrix[:,betapos:(betapos + k[j].shape[1])] * k[j])
            if betascore > alphascore:
                alphascore = betascore
                alphapos = betapos
    
        except ValueError:
            pass
    secondscoredic[j] = alphascore

print sum(topscoredic.values())
print sum(secondscoredic.values())    
print sum(worstscoredic.values()) 


  



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
    gene = str.replace(j,"Homosapiens_CISBP/pwmsUnif/","")    
    gene = str.replace(gene,".txt","")
    genelist.append(gene)
    faildic[gene] = 0 
    testno += 1
    #print testno, gene

# print sum(topscoredic.values())

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
#c = matrixmaker(startsequence)

c = randomseq(n,"numpy")


## generate a list of all seeded promoter files inthe output directory to determine which number to add for this output. 
promoterlist = []
for i in os.listdir("promoter_files/seeded/"+str(n)+"bp/"):
    if i.endswith(".txt"):
        promoterlist.append("promoter_files/seeded/"+str(n)+"bp/"+i)

pversions = [] # lsit of promoter version numbers from the promoterlist
for i in promoterlist:
    x = i
    x = x.replace("promoter_files/seeded/"+str(n)+"bp/montecarlosequence"+str(n)+"bp_", "")
    x = x.replace(".txt","")
    pversions.append(int(x))

f  = max(pversions) + 1 # current version for this promoter run.

with open("promoter_files/seeded/"+str(n)+"bp/montecarlosequence"+str(n)+"bp_" + str(f)+ ".txt", "w") as mcseq:
    mcseq.write("")

bestalign = 0 
notlist = genelist




while loop < 10000000:
    
    posdic = {}
    scoredic = {}
    seqmatrix = c
    corecount = 0
    align = 0
    keylist = []
    
    if REcheck(seqmatrix, letters, pgl4RE) == 0:
        pass
    else:
        for j in k:    
            alphapos = 0
            alphascore = -1000    
            gene = str.replace(j,"Homosapiens_CISBP/pwmsUnif/","")    
            gene = str.replace(gene,".txt","")
            genelist.append(gene)

            for i in xrange(n):
                try:
                    betapos = i 
                    betascore = np.sum(seqmatrix[:,betapos:(betapos + k[j].shape[1])] * k[j])
                    if betascore > alphascore:
                        alphascore = betascore
                        alphapos = betapos
                except ValueError:
                    pass
            


            if alphascore < float(scorethrdict[j]):
                #loop += 1
                #c = promotergenerator(n,corelist,corespec, "numpy")
                #print corecount, j, alphascore, float(scorethrdict[j]), "N"
                notalign = list(set(genelist).difference(keylist))
                faildic[gene] = faildic[gene]+1
                if corecount > bestalign:
                    bestalign = corecount 
                    notlist = notalign 
                #print loop, gene, alphapos, alphascore
                # print gene, corefinal[gene]
                # print matrixconverter(seqmatrix,letters)[alphapos:alphapos+k[j].shape[1]:]
                # if loop % 50 == 0:
                #     print loop, omegascore, testalign   
                break

            else:
                corecount += 1
                scoredic[j] = alphascore
                keylist.append(gene)
                
                #print corecount, j, alphascore, float(scorethrdict[j]), "Y"


            # for m in corefinal[gene]:
            #     print gene, m
            #     print matrixconverter(seqmatrix,letters)[alphapos:alphapos+k[j].shape[1]:]
            #     fitcount = 0
            #     if m in matrixconverter(seqmatrix,letters)[alphapos:alphapos+k[j].shape[1]:]: # test if the best alignment has the core motif. 
            #         print gene, m
            #         print matrixconverter(seqmatrix,letters)[alphapos:alphapos+k[j].shape[1]:]
            #     #if m in matrixconverter(seqmatrix,letters):
            #         fitcount += 1
            #     else: 
            #         print gene, m
            #         print matrixconverter(seqmatrix,letters)[alphapos:alphapos+k[j].shape[1]:] , "NOT"


            # if fitcount == 0:
            #     notalign = list(set(genelist).difference(keylist))
            #     faildic[gene] = faildic[gene]+1
            #     if corecount > bestalign:
            #         bestalign = corecount 
            #         notlist = notalign 
            #     break

            # else:
            #     corecount += 1
            #     scoredic[j] = alphascore
            #     keylist.append(gene)
                
            

        
      
      
            


        if sum(scoredic.values()) > omegascore and corecount == 76:
            
            bestalign = len(k)
            notlist = []
            omegascore = sum(scoredic.values())
            seqfinal = matrixconverter(seqmatrix, letters)
            omegascoredic = scoredic
            omegaposdic = posdic
            omegaalign = corecount
            with open("promoter_files/seeded/"+str(n)+"bp/montecarlosequence"+str(n)+"bp_" + str(f)+ ".txt", "w") as mcseq:
                mcseq.write(seqfinal)


    loop += 1
    
    c = promotergenerator(n,corelist,corespec, "numpy")


    #output = "%s,%s\n" %(loop, omegascore) 
    #with open("montecarlolog.txt", "a") as matlablog:
    #    matlablog.write(output)


    
    if loop % 50 == 0:
        print loop, omegascore, bestalign, n, notlist, seqfinal
   
   



##### Problem: Regenerateseq is not working properly. and ideally, the regenerate LISTS should be called by promoter functions, not the regenerate sequence itself. Otherwise, the computer needs to remake the list at every loop/. 
### After some testing, the regenerateseq is working. but it doesnt give good promoters. My hypothesis is that using core sequences is too rigid. Never has a chance to be aligned. 
## Generate a score distribution for each PWM. However, do NOT generate the score from every possible sequence for that PWM, but only all posible sequences using the core motif and random residues at other positions. 



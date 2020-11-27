# -*- coding: utf-8 -*-

import copy
import math
import sys
import os

def Kmers(k):
    kmers = [] 
    oligos = ['A','G','C','T']
    if k == 1:
        kmers = oligos
        return kmers  
    k_1mers = Kmers(k-1) 
    for oligo in oligos:
        for k_1mer in k_1mers:
            kmer = k_1mer + oligo
            kmers.append(kmer)
    return kmers

#built overall position weight matrix（allsiteNum=42-k）
def PWM(sequenceFile,kmers,allsiteNum,p0):
    pwm = []
    k = len(kmers[0])
    realCounts = []
    for i in range(allsiteNum):
        name = {}
        for kmer in kmers:
            name[kmer] = 0
        realCounts.append(name)
    pwm = copy.deepcopy(realCounts)
    totalCounts = []
    psetotalCounts = []
    ##############count the number of kmers at the corresponding allsites
    fobjr = open(sequenceFile,'r')
    for eachline in fobjr:
        if eachline[0] == '>':
            pass
        else:
            eachline = eachline.strip()
            length = len(eachline)
            for i in range(length+1-k):
                kmer = eachline[i:i+k]
                realCounts[i][kmer] += 1                          
    fobjr.close()
    for i in range(allsiteNum):
        totalCounts.append(sum(realCounts[i].values()))    ###calculate sum of all values of dictionary
        #print("the total number of real counts at the %d site is %d" % (i,sum(realCounts[i].values())))
        psetotalCounts.append(math.sqrt(totalCounts[i]))   ###math.sqrt() calculate square root
        #print("psetotalCounts at the %d site is %d" % (i,psetotalCounts[i]))
    for i in range(allsiteNum):
        for kmer in kmers:
            pwm[i][kmer] = (realCounts[i][kmer] + p0 * psetotalCounts[i])/(totalCounts[i]+ psetotalCounts[i])
    #print('the size of pwm is %d' % (len(pwm)))
    return pwm

def PCSF(sequences,kmers,p0,pwm,sites):
    pcsf = []
    k = len(kmers[0])
    for eachline in sequences:
        if eachline[0] == '>':
            pass
        else:
            eachline = eachline.strip()
            pcsfI = []
            length = len(eachline)
            for i in range(length+1-k):
                if i in sites:
                    kmer = eachline[i:i+k]
                    pcsfI.append(math.log(pwm[i][kmer]/p0,math.e))
                else:
                    pass             
            pcsf.append(pcsfI)
    #print('the size of pcsf is %d' % (len(pcsf)))
    return pcsf

def writePcsfFea(pcsf,num,pcsfFeaFile):    
    fobjw = open(pcsfFeaFile,'w')
    for i in range(len(pcsf)):     
        for j in range(len(pcsf[i])):
            temp = num + j        
            fobjw.write('%d:%f\t' % (temp,pcsf[i][j]))
        fobjw.write('\n')		
    fobjw.close()
        
def getPCSF(sequences,pcsfFeaFile):
    k = 4
    kmers = Kmers(k)
    allsiteNum = 42-k+1 
    p0 = 1/4**k
    sites = [10,11,12,13,14,15,16,17,18,19,20,21,22,23] #conservation sites by conservation.py and sorts of analysis.
    pwmP = PWM('inputFile.fasta',kmers,allsiteNum,p0) 
    pcsf = PCSF(sequences,kmers,p0,pwmP,sites) 
    writePcsfFea(pcsf,1,pcsfFeaFile)

'''
sequences = ['AGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAATAAAGAGAGA','CGGGCCTATAAGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAA','TATGTAACATAATGCGACCAATAATCGTAATGAATATGAGAAGTGTGATATTATAACATTTCATGACTACTGCAAGACTAA','TCGCACGGGTGGATAAGCGTTTACAGTTTTCGCAAGCTCGTAAAAGCAGTACAGTGCACCGTAAGAAAATTACAAGTATAC']

getPCSF(sequences,'pcsfFeaFile.txt')'''

if __name__ == '__main__':
    #sequences = ['CGCTCTATCCTGGGTTTTTGGCTGTGCCAAAAGGGAATAATGAAAAACAATAGCATCTTTGTGAAGTTTGTATTATAATAA','GTTTCCCTTATTTTTTGATAAAAGGCTTCCGAAGAAACGTAACTGTGGTATGATGTATGGAAGATAGCTAGGAACGGATTT','GGTCCGTTTTTTTGTAAAAAAAGGATTGACTTTGTGAGTCAAAGTATTTATTGTATTAAGTGTACTAATTGAAGTAATACA','GTATATATTATTTTTCATGTTTAGACAATTTTCGTCAAATTATTTGATATACTTAGGGGTGAAAGCCGCGCGTATTGTAAG']
    #read file and save list
    inputFile = open(sys.argv[1])
    sequences = []

    for eachline in inputFile:
        eachline = eachline.strip('\n')
        if eachline[0] == '>':
                pass
        else:
               sequences.append(eachline)
               
    inputFile.close()

    ##define outputFile name
    inputFilename = os.path.basename(sys.argv[1])
    #outputFilename = os.path.basename(sys.argv[2])
    
    getPCSF(sequences, 'PcsfFea_'+inputFilename)

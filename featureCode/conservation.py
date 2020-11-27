# -*- coding: utf-8 -*-

def conservation(allSequenceFile,kmers):  
    mk = []                              #the overall conservation of position i
    mks = []                             #conservation value of every Kmer at position i
    probabilities = []                   #frequency value of every Kmer at position i
    counts = []                          #count number of every Kmer at position i
    k = len(kmers[0])
    print("k = %d"%(k))
    pe = 1/(4**k)                    #background frequency
    sampleNum = 0
    for i in range(42-k+1):
        siteI = {}
        for kmer in kmers:
            siteI[kmer] = 0
        counts.append(siteI)
    print('the size of sites is %d' % (len(counts)))    
    probabilities = copy.deepcopy(counts)
    mks = copy.deepcopy(counts)
    fobjr = open(allSequenceFile,'r')
    for eachline in fobjr:
        if eachline[0] == '>':
            pass
        else:
            sampleNum += 1
            eachline = eachline.strip()
            #print('the length of eachline is %d' % (len(eachline)))
            length = len(eachline)
            for i in range(length+1-k):
                kmer = eachline[i:i+k] 
                counts[i][kmer] += 1

    fobjr.close() 
    print('the number of samples is %d' % (sampleNum))
    for i in range(42-k):
        for kmer in kmers:
            #probabilities[i][kmer] = counts[i][kmer]/sampleNum
            probabilities[i][kmer] = counts[i][kmer]/sum(counts[i].values())
            mks[i][kmer] = (probabilities[i][kmer]-pe)**2/pe
        mk.append(sum(mks[i].values()))
        """
    for i in range(301-k):
        for kmer in kmers:
            mks[i][kmer] = (probabilities[i][kmer]-pe)**2/pe
        mk.append(sum(mks[i].values()))
        """
    return mk

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

def writeConservation(filename,mk):
    fobjw = open(filename,'w')
    for i in range(len(mk)):
        fobjw.write("%d\t%f\n" %(i+1,mk[i]))
    fobjw.close()


import copy
import sys
sequenceFile = sys.argv[1] 
outputFile = sys.argv[2]

if __name__=='__main__':
    
    for i in range(4,5):
        kmers = Kmers(i)
        mk = conservation(sequenceFile,kmers)
        #filename = '%dmers-promoter.txt' % (i)
        writeConservation(outputFile,mk)

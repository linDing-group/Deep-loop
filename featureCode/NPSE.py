# -*- coding: utf-8 -*-
'''
@author: dfy

'''

def CalculateFrequency(input_seq_file, gap):
	line_num = 0#calculate how many sequences
	frequency = []
	seq_name = []
	all_kmer = fourth_file(2)
	f = open(input_seq_file,'r')

	lines = f.readlines()
	for line in lines:
		if line[0] == '>':
			seq_name.append(line.strip('\n').replace('>',''))
		else:
			line = line.strip()
			seqLen = len(line)
			kmer = []

			for i in range(seqLen-gap-1):
				kmer.append('%s%s'%(line[i],line[i+gap+1]))
			#print(kmer)

			kmer_len = len(kmer)	
			each_fre = []
			for each in all_kmer:
				fre = kmer.count(each)/kmer_len
				each_fre.append(fre)
			frequency.append(each_fre)
			line_num += 1
	len_kmer = len(all_kmer)
	return frequency,all_kmer,line_num,len_kmer,seq_name


def fourth_file(k):
	a = []
	b = ['A','G','C','T']
	if k == 2:
		for i in b:
			for j in b:
				a.append(i+j)
		return a

outfile = '.txt' #outputFile
f = '.fasta' #inputfile
gap = 4 #here is the parameter of max gap
frequency,all_kmer,line_num,len_kmer,seq_name = CalculateFrequency(f,gap)
#print(frequency)

g = open(outfile,'w')
for num in range(line_num):
	line_fea = ''
	for num2 in range(16):
		if num2 == 15:
			line_fea += '%d:%.6f\n'%(num2+1,frequency[num][num2])
		else:
			line_fea += '%d:%.6f\t'%(num2+1,frequency[num][num2])
	line_fea = '1' +'\t'+line_fea
	g.write(line_fea)



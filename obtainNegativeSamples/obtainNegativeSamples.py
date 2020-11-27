
import copy
import math


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
        totalCounts.append(sum(realCounts[i].values()))
        #print("the total number of real counts at the %d site is %d" % (i,sum(realCounts[i].values())))
        psetotalCounts.append(math.sqrt(totalCounts[i]))
        #print("psetotalCounts at the %d site is %d" % (i,psetotalCounts[i]))
    for i in range(allsiteNum):
        for kmer in kmers:
            pwm[i][kmer] = (realCounts[i][kmer] + p0 * psetotalCounts[i])/(totalCounts[i]+ psetotalCounts[i])
    #print('the size of pwm is %d' % (len(pwm)))
    return pwm
def threValue(filename):
	k = 1
	sites = []
	kmers = Kmers(k)
	allsiteNum = 42-k+1 
	p0 = 1/4**k
	for i in range(allsiteNum+1):
	    sites.append(i)
	sequenceFile = filename
	pwmP = PWM(sequenceFile,kmers,allsiteNum,p0)
	#print(len(pwmP))

	f1 = open(filename,'r')
	#g = open('neg_R.txt','w')
	p0 = 1/4
	window = 42
	all_v = []
	for m in f1:
	    m = m.strip('\n')
	    if '>' in m:
	        #lable = m.strip('\n')
	        chrL = m.split(':')[0].replace('>','')
	        s = int(m.split('-')[0].split(':')[-1])
	    else:
	        value_pos = []
	        subSeq_pos = []
	        mm = m.upper()
	        for n in range(0,len(mm)-window+1):
	            value_each_pos = 0
	            o = 0
	            subSeq_pos.append(mm[n:n+window])
	            for each in mm[n:n+window]:
	                #value_each_pos += PWM[o][each]
	                value_each_pos += math.log(pwmP[o][each]/p0,math.e)
	                o += 1
	            value_pos.append(value_each_pos)
	        Ind_pos = value_pos.index(max(value_pos))
	        #print(max(value_pos))
	        all_v.append(max(value_pos))
	return(min(all_v))

def getRangeSeq(fname,g1name, g2name):
	f = open(fname,'r')
	g1 = open(g1name,'w')
	g2 = open(g2name,'w')
	for i in f:
		if '>' in i:
			chro1 = i.split('\t')[0].split(':')[0].split('...')[-1]
			# chro2 = i.strip('\n').split('\t')[-1].split(':')[0].split('...')[-1]
			start = i.split('\t')[0].split(':')[-1].split('-')[0]#first sequence start
			end = i.strip('\n').split('\t')[-1].split(':')[-1].split('-')[0]#second sequencestart
			# dis = abs(int(end)-int(start))

			e1 = int(start)-1100
			s1 = e1-100
			e2 = int(end)-1100
			s2 = e2-100
			g1.write('%s\t%d\t%d\n'%(chro1,s1,e1))
			g2.write('%s\t%d\t%d\n'%(chro1,s2,e2))

			s11 = int(start)+1100
			e11 = s11+100
			s22 = int(end)+1100
			e22 = s22+100
			g1.write('%s\t%d\t%d\n'%(chro1,s11,e11))
			g2.write('%s\t%d\t%d\n'%(chro1,s22,e22))

def final_out_file(a_file, b_file, f_filename, g_filename):
	a = threValue(a_file)
	b = threValue(b_file)
	f = open(f_filename,'r')
	g = open(g_filename,'w')
	ff = f.readlines()
	for i in range(len(ff)):
		if '>' in ff[i]:
			ii = ff[i].strip('\n').split('\t')
			if float(ii[1]) > a and float(ii[-1]) > b:
				g.write(ff[i]+ff[i+1])

def obtain_pwm_neg(seqfile, f1name, gname):
    k = 1
    sites = []
    kmers = Kmers(k)
    allsiteNum = 42-k+1 
    p0 = 1/4**k
    for i in range(allsiteNum+1):
        sites.append(i)
    sequenceFile = seqfile
    pwmP = PWM(sequenceFile,kmers,allsiteNum,p0)
    #print(pwmP)

    f1 = open(f1name,'r')
    g = open(gname,'w')
    p0 = 1/4
    window = 42
    all_v = []
    for m in f1:
        m = m.strip('\n')
        if '>' in m:
            #lable = m.strip('\n')
            chrL = m.split(':')[0].replace('>','')
            s = int(m.split('-')[0].split(':')[-1])
        else:
            value_pos = []
            subSeq_pos = []
            mm = m.upper()
            for n in range(0,len(mm)-window+1):
                value_each_pos = 0
                o = 0
                subSeq_pos.append(mm[n:n+window])
                for each in mm[n:n+window]:
                    #value_each_pos += PWM[o][each]
                    value_each_pos += math.log(pwmP[o][each]/p0,math.e)
                    o += 1
                value_pos.append(value_each_pos)
            Ind_pos = value_pos.index(max(value_pos))
            #print(max(value_pos))
            all_v.append(max(value_pos))
            # print(subSeq_pos[Ind_pos])
            g.write('>%s:%d-%d\t%f\n%s\n'%(chrL,s+Ind_pos,s+Ind_pos+window-1,max(value_pos),subSeq_pos[Ind_pos]))
def getFinalFa(f1name,f2name,gname,f1_posname,f2_posname,finaloutfile):
	f1 = open(f1name,'r')
	f2 = open(f2name,'r')
	g = open(gname,'w')
	ff1 = f1.readlines()
	ff2 = f2.readlines()
	print(len(ff1),len(ff2))
	for i in range(0,len(ff1)-1,2):
		g.write(ff1[i].strip('\n')+'\t'+ff2[i]+ff1[i+1].strip('\n')+'\t\t'+ff2[i+1])

	final_out_file(f1_posname, f2_posname, gname, finaloutfile)

def split(fname,g1name,g2name):
	f = open(fname,'r')
	g1 = open(g1name,'w')
	g2 = open(g2name,'w')
	for i in f:
		if 'Orien' in i:
			pass
		else:
			ii = i.split('\t')
			g1.write(ii[0]+'\n')
			g2.write(ii[1])
	f.close()
	g1.close()
	g2.close()


#step1: According to the coordinate of positive sample to obtain negative sample's, 
# *-neg.txt was the coordinate of negative sample,
# *-neg.fa was obtained by BedTools
inputname = 'FR'
split(inputname+'.txt',inputname+'-F.txt',inputname+'-R.txt')
getRangeSeq(inputname+'.txt',inputname+'-F-neg.txt', inputname+'-R-neg.txt',)

#step2:According to the PWM of ositive sample, the 42bp sequence of negative sample with the largest score was taken out and wrote in file *neg-final.fa

obtain_pwm_neg(inputname+'-R.txt', inputname+'-R-neg.fa', inputname+'-R-neg-final.fa')
obtain_pwm_neg(inputname+'-F.txt', inputname+'-F-neg.fa', inputname+'-F-neg-final.fa')

#step3:'FR-neg-final_out.fa' is finally obtained negative sample sequences.
getFinalFa(inputname+'-F-neg-final.fa',inputname+'-R-neg-final.fa',inputname+'-neg-final.fa',inputname+'-F.txt',inputname+'-R.txt',inputname+'-neg-final_out.fa')


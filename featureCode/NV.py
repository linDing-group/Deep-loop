# coding=utf-8

import os
import multiprocessing

###convert every A, C, F, and T to number 0, 1, 2, and 3###
def tran_digital(ch):
    if ch == 'A':
        return 0
    elif ch == 'C':
        return 1
    elif ch == 'G':
        return 2
    else:
        return 3


##convert every sequence to decimal number. Foe example, if window size is 3, AAAAAA-TTTTTT could be represent 0-4095.###
def cal_label(feature):
    index = 0
    for i in range(0, len(feature)):
        index = 4 * index + feature[i]
    return index

win_size = 1 #window 
filename = 'nature_R2_' + str(win_size) 
out = open(filename + ".txt", "w") #outputFile

###read positive samples
f1 = open('.txt', 'r')
for line in f1:
    if line[0] == '>':  # pass annotation line
        continue
    line = line.upper()  # convert every lower nucleotide to upper nucleotide
    line = line.strip('\n')  
    number = [0] * len(line)

    # convert every A, C, F, and T to number 0, 1, 2, and 3
    for i in range(0, len(line)):
        number[i] = tran_digital(line[i])

    # calculate frequency
    nk_frequence = [0] * (4 ** win_size)
    tk_position = [0] * (4 ** win_size)
    u = [0] * (4 ** win_size)
    for i in range(0, len(number) - win_size + 1):
        id = cal_label(number[i:i + win_size])
        nk_frequence[id] += 1
        tk_position[id] +=i+1
    for i in range(0,4**win_size):
        if nk_frequence[i] == 0 :
           u[i] = 0
        else:
           u[i] = tk_position[i]/nk_frequence[i]

	
    dk_secondery = [0] * (4 ** win_size)
    for i in range(0, len(number) - win_size + 1):
        id = cal_label(number[i:i + win_size])
        dk_secondery[id] += ((i+1-u[id])**2)/(nk_frequence[id]*len(number))
		
    out.write("1\t")
    for i in range(0,len(dk_secondery)):
        out.write("%d:%.4f\t" % (1+3*i, nk_frequence[i]))
        out.write("%d:%.4f\t" % (2+3*i, u[i]))
        out.write("%d:%.4f\t" % (3+3*i, dk_secondery[i]))
    out.write('\n')
	
###read negative samples
f1 = open('.txt', 'r')
for line in f1:
    if line[0] == '>':  
        continue
    line = line.upper() 
    line = line.strip('\n') 
    number = [0] * len(line)

    for i in range(0, len(line)):
        number[i] = tran_digital(line[i])

    # calculate frequency
    nk_frequence = [0] * (4 ** win_size)
    tk_position = [0] * (4 ** win_size)
    u = [0] * (4 ** win_size)
    for i in range(0, len(number) - win_size + 1):
        id = cal_label(number[i:i + win_size])
        nk_frequence[id] += 1
        tk_position[id] +=i+1
    for i in range(0,4**win_size):
        if nk_frequence[i] == 0 :
           u[i] = 0
        else:
           u[i] = tk_position[i]/nk_frequence[i]

	
    dk_secondery = [0] * (4 ** win_size)
    for i in range(0, len(number) - win_size + 1):
        id = cal_label(number[i:i + win_size])
        dk_secondery[id] += ((i+1-u[id])**2)/(nk_frequence[id]*len(number))
		
    out.write("2\t")
    for i in range(0,len(dk_secondery)):
        out.write("%d:%.4f\t" % (1+3*i, nk_frequence[i]))
        out.write("%d:%.4f\t" % (2+3*i, u[i]))
        out.write("%d:%.4f\t" % (3+3*i, dk_secondery[i]))
    out.write('\n')
out.close()
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
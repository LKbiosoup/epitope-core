# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:52:17 2023

@author: likun
"""
'''
    This script is used to obtain the immune epitope core. 
    The modifiable lines are 134-136, with two input files and one output file path in each of these three lines
 134   file='D:/imm/file1.csv' #This line needs to specify a CSV file containing all epitopes, which can be found in the https://github.com/LKbiosoup/epitope-core named file1
 135   file2='D:/imm/file2.csv'#This is a CSV file containing 5000 epitopes, used to detect the coverage of the obtained epitope core to the epitope, which can be found in the https://github.com/LKbiosoup/epitope-core named file2
 136   filew='D:/imm/run.txt'  #This file records the running results, and scripts typically take at least a few days to complete. When the coverage tends to balance, it can interrupt the operation in advance
    #In the record file(run.txt) ,the last three lines record the latest running results. The first line of these three lines records the core obtained, and the second line records the coverage. Report: [0.0552, 0.1888, 0.2877, 0.2743, 0.194] corresponds to the proportion of detected epitope cores with 0, 1, 2, 3, or more in the epitope library   
     #In the running results, the list named setx records all acquired Epitope cores 
 '''

# Import libraries
from Bio import AlignIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
from threading import Thread
#from Bio.SubsMat import MatrixInfo 
from Bio.Cluster import kcluster
import difflib
from threading import Thread


def jb(seq1,seq2):
    try:
        matrix = substitution_matrices.load("BLOSUM62")
        alignments = pairwise2.align.localds(seq1, seq2, matrix, -5, -2)
        best_alignment = alignments[0]
        best_similarity=best_alignment.score
        return  best_similarity
    except:
        print()

def cl(seq):
    seq=seq.replace('"','')
    seq=seq.split('+')[0]
    seq=seq.strip()
    return seq

 
def treadx(seq1,lines):

    top1000=['s']*100;top1000s=[0]*100
    min_index=0
    for i in lines:
        seq2=cl(i)
        score=jb(seq1,seq2)  
        try:
            if score>top1000s[min_index]: 
                if seq2 not in seq1 and seq1 not in seq2:
                    top1000s[min_index]=score
                    top1000[min_index]=seq2
                    min_index = top1000s.index(min(top1000s))
        except:
            print('error line')
    lis= [[top1000[x],top1000s[x]] for x in range(len(top1000))]
    lis=sorted(lis, key=lambda x: x[1], reverse=True)
    return lis

def longest_common_subsequence(a, b,num):
    sm = difflib.SequenceMatcher(None, a, b)

    match = sm.find_longest_match(0, len(a), 0, len(b))

    if match.size >= num:

        lcs = a[match.a:match.a + match.size]
        return lcs

def core(seq1,lines):

    sets=treadx(seq1,lines)#
    set2=[sets[x][0] for x in range(len(sets))]
    cores=[]
    for seq2 in set2:
        if  longest_common_subsequence(seq1, seq2,4):
            cores.append(longest_common_subsequence(seq1, seq2,4))
            
    cores1=list(set(cores))
    #print(cores1) 
    x=[]
    for i in cores1:
        tag=1
        for j in cores1:
            if j in i and len(j)<len(i):
                tag=0
        if tag!=0:
            x.append(i)
    
    while tag==0:
        tag=1
        pairs = [[a,b] for i, a in enumerate(x) for b in x[i+1:]]
        for pair in pairs:
            seq1=pair[0];seq2=pair[1]
            if longest_common_subsequence(seq1, seq2,4):
                tag=0
                x.remove(seq1);x.remove(seq2)
                x.append(longest_common_subsequence(seq1, seq2,4))
                break             
    return x
def yz(sets,file):
    with open (file,'r') as f:
        lines=f.readlines()
    f.close()
    length=len(lines)
    out=[]
    for line in lines:
        line=cl(line)
        tag=0;set2=[]
        for x in sets:
            if x in line:
                if tag==0:
                    tag=tag+1
                    set2.append(x)
                else:
                    if len([set2[i] for i in range(len(set2)) if longest_common_subsequence(x,set2[i],3)==None])==len(set2):
                        #print([set2[i] for i in range(len(set2)) if longest_common_subsequence(x,set2[i],3)==None])
                        tag=tag+1
                        set2.append(x)
        out.append(tag)
    x0=round(len([x for x in out if x==0])/length,4)
    x1=round(len([x for x in out if x==1])/length,4)
    x2=round(len([x for x in out if x==2])/length,4)
    x3=round(len([x for x in out if x==3])/length,4)
    other=round(1-x0-x1-x2-x3,4)    
    return [x0,x1,x2,x3,other]
        
if __name__ == "__main__":         
    file=''                     #all epitopes
    file2=''                    # A portion of epitopes used to detect coverage
    filew='D:/imm/run.txt'       # Record of operation results
    setx=[]                     #epitope cores result

    with open (file,'r') as f:
        lines=f.readlines()
    f.close()
    tag=1
    lines1=[lines[i] for i in range(len(lines)) if i%20==0]
    lines2=[lines[i] for i in range(len(lines)) if i%20!=0]
    for line1 in lines1:
        
        line1=cl(line1)
        print(line1)
        lin=core(line1,lines2)
        setx.extend(lin)
        
        setx=list(set(setx))
        df=yz(setx,file2)
        with open(filew,'a') as fw:
            fw.write(str(setx)+'\n')
            fw.write(' report: '+str(df)+'\n')
            fw.write(str(tag)+'\n')    
        fw.close()
        tag=tag+1
   
    
    
        

        
    

    

    



    



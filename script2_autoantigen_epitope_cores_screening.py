
"""
Created on Mon Jun 19 17:06:12 2023

@author: likun
"""
'''
This script is executed to remove redundancy after obtaining the initial antigen epitope core. The following three lines need to be changed
 line 295: runfile='E:/imm/run.txt'  # Runfile is the running record file for demo1
 line 296: file='E:/imm/file5.csv'     # File3 contains all autoantigen epitopes, which can be found in the https://github.com/LKbiosoup/epitope-core 

 #
 #/
 # In the running results, list liss4 contains the core of the autoantigen epitope and the corresponding number of occurrences in the pathogen epitope. List score represents the coverage of the core over the entire pathogen epitope library
        

'''

from Bio import AlignIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
from threading import Thread
#from Bio.SubsMat import MatrixInfo 
from Bio.Cluster import kcluster
from threading import Thread
import difflib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import lfilter
from pulp import * 
import os 
import numpy as np
def longest_common_subsequence(a, b,num):
    sm = difflib.SequenceMatcher(None, a, b)
    match = sm.find_longest_match(0, len(a), 0, len(b))

    if match.size >= num:
        lcs = a[match.a:match.a + match.size]
        return lcs
def cl(seq):
    seq=seq.replace('"','')
    seq=seq.split('+')[0]
    seq=seq.strip()
    return seq

def yz(sets,file):
    with open (file,'r') as f:
        lines=f.readlines()
    f.close()
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
    x0=round(len([x for x in out if x==0])/len(lines),4)
    x1=round(len([x for x in out if x==1])/len(lines),4)
    x2=round(len([x for x in out if x==2])/len(lines),4)
    x3=round(len([x for x in out if x==3])/len(lines),4)
    other=round(1-x0-x1-x2-x3,4)    
    return [x0,x1,x2,x3,other]



def yz2(sets,file):
    with open (file,'r') as f:
        lines=f.readlines()
    f.close()
    out=[]
    for line in lines:
        line=cl(line)
        tag=0
        for x in sets:
            if x in line:
                tag+=1
        out.append(tag)
    x0=round(len([x for x in out if x==0])/len(lines),4)
    x1=round(len([x for x in out if x==1])/len(lines),4)
    x2=round(len([x for x in out if x==2])/len(lines),4)
    x3=round(len([x for x in out if x==3])/len(lines),4)
    other=round(1-x0-x1-x2-x3,4)    
    return [x0,x1,x2,x3,other]

                    
def compress(lis,file):
    
    score=[[i,0] for i in lis ]
    with open (file,'r') as f:
        lines=f.readlines()
    f.close()
    
    for core in range(len(lis)):
        for line in lines:
            if lis[core] in line:
                score[core][1]+=1
       
    sorted_lis = sorted(score, key=lambda x: x[1], reverse=True)  
    sorted_lis=sorted_lis[0:len(sorted_lis)]
    sorted_lis=[i[0] for i in sorted_lis]
    print('bingo')
    return sorted_lis      
def compress2(lis,file):
    
    score=[[i,0] for i in lis ]
    with open (file,'r') as f:
        lines=f.readlines()
    f.close()
    
    for core in range(len(lis)):
        for line in lines:
            if lis[core] in line:
                score[core][1]+=1
       
    sorted_lis = sorted(score, key=lambda x: x[1], reverse=True)  
    sorted_lis=sorted_lis[0:len(sorted_lis)]
    #sorted_lis=[i[0] for i in sorted_lis]
    print('bingo')
    return sorted_lis  
                    
def amp(sets):
    amset='ACDEFGHIKLMNPQRSTVWY'
    out=[]
    for am in amset:
        num=0
        for core in sets:
            if am in core:
                num+=1
        x=[am,num/len(sets)]
        out.append(x)
    return out
def amp3(sets,am):
    num=0
    for core in sets:
        tag=0
        for amo in am:
            
            if amo in core:
                tag=1
        if tag==1:
            num+=1
    return num/len(sets)

    

def out_ry3(sets,file,outfile):
   
    with open(file,'r') as f:
        lines=f.readlines()
    f.close()
    out=[]
    for core in sets:
        coresp=[]
        for line in lines:
            seq=line.split(',')[0]
            if core[0] in seq:
                sp=line.split(',')[1]
                sp=sp.strip()
                if 'SARS' in sp:
                    sp='SARS'
                elif 'Human coronavirus' in sp:
                    sp='Human coronavirus'
                if sp not in coresp and len(sp)>3:
                    coresp.append(sp)
        if len(coresp)>3:
            out.append(core)
            with open(outfile,'a') as f:
                f.write(str(core)+'\n')
                f.write(str(coresp)+'\n')
            f.close()
    return out 
def concatenate(s1, s2):
    com=longest_common_subsequence(s1,s2,len(s1)-1)
    ind1=s1.index(com)
    ind2=s2.index(com)
    if ind1 == ind2:
        return False
    elif ind1==0 and ind2==1:
        return s2+s1[-1]
    elif ind1==1 and ind2==0:
        return  s1+s2[-1]
def getcon(sets,file):
    out=[]
    tag=0
    for i in sets:
        lf=[x for x in sets if len(x)==len(i) and i!=x and  longest_common_subsequence(i,x,len(i)-1)]
            
        if len(lf)>0:
            for l in lf:
                if concatenate(i,l):
                    lis1=[];lis2=[]
                    with open(file,'r') as f:
                        lines=f.readlines()
                    f.close()
                    for line in range(len(lines)):
                        if i in lines[line]:
                            lis1.append(line)
                        if l in lines[line]:
                            lis2.append(line)
                    if max(len(lis1),len(lis2))*0.6<len([h for h in lis1 if h in lis2]):
                        out.append(concatenate(i,l))
                        print('two core belong:',concatenate(i,l))
                        '''
                        try:
                            ind1=out.index(i)
                            del out[ind1]
                            ind2=out.index(l)
                            del out[ind2]
                        except:
                            print()
                    '''
                     
        tag+=1
        #print(tag)
    out=list(set(out))
    print('bingo')
    return out
    
def gen_substr(s, n):
    result = []
    for i in range(len(s) - n + 1):
        result.append(s[i:i+n])
    return result 
def have_same_chars(A, B):
    count = 0
    for a, b in zip(A, B):
        if a == b:
            count += 1

    return count
def maxpp(A,B,num):
    lisa=gen_substr(A,num)
    lisb=gen_substr(B,num)
    maxzf=0
    for suba in lisa:
        for subb in lisb:
            zf=have_same_chars(suba,subb)
            if maxzf>zf:
                maxzf=zf
    return maxzf
            
        
def out_ry2(sets):
    setscopy=sets[:]
    for core in sets: 
        if len(core[0])==4:
           
            sons=[x for x in sets if core[0] in x[0] and x != core]
            for son in sons:
                try:
                    if 0.15*core[1]>son[1]:
                        ind=setscopy.index(son)
                        del setscopy[ind]
                    elif 0.7*core[1]<son[1]:
                        ind=setscopy.index(core)
                        del setscopy[ind]
                except:
                        print('')
                    
        else:
            
            sons=[x for x in sets if core[0] in x[0] and x != core]
            for son in sons:
                try:
                    if 0.15*core[1]>son[1]:
                        ind=setscopy.index(son)
                        del setscopy[ind]
                    elif 0.7*core[1]<son[1]:
                        ind=setscopy.index(core)
                        del setscopy[ind]
                except:
                                print('')
            parents=[x for x in sets if x[0] in core[0] and x != core]
            for parent in parents:
                try:
                    if core[1]>parent[1]*0.7:
                        ind=setscopy.index(parent)
                        del setscopy[ind]
                    elif 0.15*parent[1]>core[1]:
                        ind=setscopy.index(core)
                        del setscopy[ind]
                except:
                    print('')
    return setscopy           
                                      
tag=1
runfile='E:/imm/run.txt'
file='E:/imm/file5.csv'

with open(runfile,'r') as f:
        lines=f.readlines()
f.close()
lines=lines[-4:]
for line in lines[-3:]:
        if '[' in line[0:3]:
            lis=line.split("', '")
            lis[0]=lis[0].replace('[','')
            lis[0]=lis[0].replace("'",'')
            lis[-1]=lis[-1].replace(']','')
            lis[-1]=lis[-1].replace("'",'')
lis=[i.strip() for i in lis]
#print(lis)
lis=compress2(lis,file)
lis5=out_ry2(lis)


lis5=[i[0] for i in lis5]

lis6=getcon(lis5,file)
lis7=getcon(lis6,file)
lis8=getcon(lis7,file)
liss=[]
for i in lis7:
    if not any(i in f for f in lis8):
        liss.append(i)
liss.extend(lis8)
liss2=[]
for i in lis6:
    if  not any(i in f for f in liss):
        liss2.append(i)
liss2.extend(liss)
liss3=[]
for i in lis5:
    if  not any(i in f for f in liss2):
        liss3.append(i)
liss3.extend(liss2)
liss4=compress2(liss3,file)
liss5=[i[0] for i in liss4 if i[1]>=5]
score=yz2(liss5,file)



      



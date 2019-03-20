#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

"""
This script takes in a allValidPairs file from hic-pro, and count interactions between specific genome loci/contig (target) and 100kb/50kb window of all chromosomes.
The count in the output files is normalized by the length of target loci/contig (per 100kb) and normalized by the length of the region it interactes with(per 100kb).
example usage: python3 count_interaction_by_window.py S_allValidPair chomosome.size chromosome_window, "2R" out.summary out.txt
"""

import sys

infile=sys.argv[1]
sizefile=sys.argv[2]
windowfile=sys.argv[3]
target=sys.argv[4]
sumfile=sys.argv[5]
outfile=sys.argv[6]


def check_range(i, start, end):
    """
    Function: check_range

    This function checks if i is in range(start, end)

    Args:
        i        
        start 
        end 

    Returns:
        result: True if i is in range, False if i is not in range
    """
    if i>=start and i<=end:
        result=True
    else:
        result=False
    return result


def get_interaction(infile, sizefile, windowfile, target, position=False):
    """
    Function: get_interaction

    This function returns a dictionary of interaction data from allValidPairs file from hic-pro.

    Args:
        infile: allValidPairs file from hic-pro
        
        sizefile: a file with contig size 
        
        windowfile: a file with 100kb window of chromosomes
        
        target: contig name you are interested in

        position: interested position range in the contig, default is false

    Returns:
        count: dict of interaction count
        summary: dict of interaction info.
        right: read name of target interaction
        other: read name of non-target interaction
    """

            
    windows={}
    with open(windowfile,'r') as wfile:
        lines=wfile.readlines()
        for line in lines:
            line=line.strip('\n')
            info=line.split('\t',)
            if info[0] not in windows.keys():
                windows[info[0]]=[]
            windows[info[0]].append(info[1]+'-'+info[2])

    #get contig sizes
    lens={}
    with open(sizefile,'r') as sfile:
        lines=sfile.readlines()
        for line in lines:
            info=line.split('\t',)
            lens[info[0]]=int(info[1])           
            
    #get the size of target loci/contig
    if position:
        start=int(position.split('-',)[0])
        end=int(position.split('-',)[1])
        size=end-start+1
    else:
        size=lens[target]
        
    summary={}
    right=[]
    other=[]    
    with open (infile, 'r') as ifile:
        lines=ifile.readlines()
        for line in lines:
            info=line.split('\t',)
            name=info[0]    
            chr1=info[1]
            pos1=int(info[2])
            chr2=info[4]
            pos2=int(info[5])

            if chr1 not in summary.keys():
                summary[chr1]={}
            if pos1 not in summary[chr1].keys():
                summary[chr1][pos1]=[]

            if chr2 not in summary.keys():
                summary[chr2]={}       
            if pos2 not in summary[chr2].keys():
                summary[chr2][pos2]=[]

            #if position set, check if pos1/pos2 in range; otherwise set pos1/pos2 to true         
            if position:
                result_pos1=check_range(pos1, start, end)
                result_pos2=check_range(pos2, start, end)
            else:
                result_pos1=True
                result_pos2=True
            
            # if both reads fall in the target contig, check if either one pos in range
            if chr1==chr2==target:
                if result_pos1:
                    summary[chr2][pos2].append((chr1,pos1))
                    right.append(name)
                elif result_pos2:
                    summary[chr1][pos1].append((chr2,pos2))
                    right.append(name)
                else:
                    other.append(name)
            # if only read1 fall in the target contig, check if pos1 in range
            elif chr1==target:
                if result_pos1:
                    summary[chr2][pos2].append((chr1,pos1))
                    right.append(name)
                else:
                    other.append(name)
            # if only read2 fall in the target contig, check if pos2 in range
            elif chr2==target:
                if result_pos2:
                    summary[chr1][pos1].append((chr2,pos2))
                    right.append(name)
                else:
                    other.append(name)
            else:
                other.append(name)
    
    # divide by windows
    middle={}
    count={}
    output={} 
    
    for key in summary.keys():
        if summary[key]!={}:
            middle[key]={}
            for win in windows[key]:
                middle[key][win]=[]
                s=int(win.split('-',)[0])
                e=int(win.split('-',)[1])
                for i in summary[key].keys():
                    # s<=i<e
                    if i in range(s,e):
                        middle[key][win].append(summary[key][i])
    
    for key in middle.keys():
        if middle[key]!={}:
            output[key]={}
            count[key]={}
            for win in middle[key].keys():
                if middle[key][win]!=[]:
                    output[key][win]=middle[key][win]
                    j=0
                    for i in middle[key][win]:
                        j+=len(i)
                    length=int(win.split('-',)[1])-int(win.split('-',)[0])
                    #raw count, and normalize by both length
                    count[key][win]=(j,100000*100000*j/size/length)
                                                                  
    return count, output, right, other



contig=target
count=get_interaction(infile, sizefile, windowfile, contig)[0]
output=get_interaction(infile, sizefile, windowfile, contig)[1]
  
 
with open(outfile,'w') as ofile:
    ofile.write("#Interaction with " + contig + "\n")
    ofile.write("%-20s\t%-20s\t%-35s\t%-20s\n" %("#Contig_1","#Position_1","#Contig_2","#Position_2"))
    for key in output.keys():
        for i in output[key].keys():
            for item in output[key][i]:
                for j in range(0,len(item)):
                    ofile.write("%-20s\t%-20s\t%-35s\t%-20s\n" %(item[j][0],item[j][1],key,i))

with open(sumfile,'w') as ofile:
    ofile.write("#Interaction with " + contig + "\n")
    ofile.write("%-20s\t%-20s\t%-20s\t%-20s\n" %("#Contig","#window","#Interactions raw count","#Interactions per 100/50kb window"))
    for key in count.keys():
        for win in count[key].keys():
            ofile.write("%-20s\t%-20s\t%-20s\t%-20s\n" %(key,win,str(count[key][win][0]),str(count[key][win][1])))


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script parse the output from Raxml, and generate input file for plotting in R.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

import re

infile="RAxML_bipartitionsBranchLabels.consensus_1_bigelements1000_with_outgroup_Jockey-3_from_Yak_and_Sec.automre"
outfile="G2_and_Jockey-3_for_R.txt"
#the list of full-length elements
full_file="full-length-G2-Jockey-3_list.txt"
#centromere contigs
cens=["3R_5","Contig119","Contig79","Y_Contig26",'tig00057289']

#parse the order of the tips
order=[]
with open(infile,'r') as ifile:
    lines=ifile.readlines()
    for line in lines:
        #split by both ( and :
        info=re.split('[(,:,)]', line)
        for item in info:
            if item.startswith('G') or item.startswith('J'):
                order.append(item)

full_list=[]
with open(full_file,'r') as ffile:
    lines=ffile.readlines()
    for line in lines:
        line=line.strip('\n')
        new=line.replace('-','_')
        full_list.append(new)

"""
pch=21 round outline -- non-cen&truncated
pch=16,round filled -- non-cen&fullength
pch=22 square outline -- cen&truncated
pch=15, square filled -- cen&fullength

"""

i=0
with open(outfile,'w') as ofile:
    
    ofile.write("%s\t%s\t%s\t%s\t%s\n" %('order','tip','color_code','name','category'))
    
    for item in order:

        #color code
        if item.startswith('J'):
            code='2'
        elif item.startswith('G'):
            code='4'
        else:
            code='0'
        
        #name three outgroup
        if item.startswith('Jockey_3_Dyak'):
            name='yak'
            code='0'
            category=''
#            fill=''

        elif item.startswith('Jockey_3_Dsec'):
            name=''
            code='3'
            category='16'
#            fill=''

        elif item == 'Jockey_3_Dsim':
            name=''
            code='5'
            category='16'
#            fill=''

        #others
        else:
            name='' 
            
            contig=item.split('.',)[1]
            # cen & full-length           
            if contig in cens and item.replace('__reversed_','') in full_list:
                category='15'               
#                fill=''
#                code=o
            # cen & truncated
            elif contig in cens:
                category='22'
                fill='white'

            # non-cen & full-length
            elif item.replace('__reversed_','') in full_list:
                category='16'
#                fill=''
#                code=o
            # non-cen & truncated
            else:
                category='21'
#                fill='white'
            
            #full length or not
            if item.replace('__reversed_','') in full_list:
                length='full_length'
            else:
                length=''            
        
        #order
        i+=1
                
        ofile.write("%s\t%s\t%s\t%s\t%s\n" %(str(i),item,code,name,category))
               

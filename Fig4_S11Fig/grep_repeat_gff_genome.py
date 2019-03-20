#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script takes in the assembly and annotation gff file, extract variant sequences based on the annotations.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

genome="dmel_scaffold2_plus0310.fasta"
gff="dmel_scaffold2_plus0310.repeat.gff3"
cens=["3R_5","Contig119","Contig79","Y_Contig26",'tig00057289']

           
TE='G2_DM'


test=set()
annotation={}
with open(gff,'r') as gfile:
    lines=gfile.readlines()
    for line in lines:
        # change == or startswith base on situation
        if not line.startswith('#'):
            info=line.split('\t',)
            # info[8] is the annotation attribute of repeats
            if info[8].split(' ',)[0].startswith(TE):
                test.add(info[8].split(' ',)[0])
                if info[0] not in annotation.keys():
                    annotation[info[0]]=[]
                length=int(info[4])-int(info[3])
                annotation[info[0]].append((info[3]+'.'+info[4]+'.'+str(length)))

# print repeat names you grep out from the gff file
print ("Repeat names\n"+test)

TE=TE.strip('_DM').strip('_DSim')
out='TE+'_variance.fasta'

#get sequence from genome assembly
seq={}
with open (genome,'r') as seqfile:
    lines=seqfile.readlines()
    for line in lines:
        line=line.strip('\n')
        if line.startswith('>'):
            i=lines.index(line+'\n')
            seq[line.strip('>')]=lines[i+1]

#output repeat variants sequences            
with open(out,'w') as ofile:
    for contig in annotation.keys():
        for item in annotation[contig]:
            start=int(item.split('.',)[0])
            end=int(item.split('.',)[1])
            length=item.split('.',)[2]
            monomer=seq[contig][(start-1):end]
            ofile.write(">%s.%s.%s-%s.%s\n%s\n" %(TE,contig,str(start),str(end),length,monomer))

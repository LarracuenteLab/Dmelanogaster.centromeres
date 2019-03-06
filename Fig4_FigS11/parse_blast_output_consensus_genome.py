#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script pulls out variant sequences based on blast result from blasting consensus to genome.

@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

genome='dere-all-chromosome-r1.3.fasta'
blast='dere.blast'
out='dere_variants.fasta'

#get genome assembly sequence
seq={}
with open (genome,'r') as seqfile:
    lines=seqfile.readlines()
    for line in lines:
        line=line.strip('\n')
        if line.startswith('>'):
            i=lines.index(line+'\n')
            header=line.split(' ',)[0].strip('>')
            seq[header]=lines[i+1]

#output repeat variants sequences            
with open(out,'w') as ofile:
    with open(blast,'r') as ifile:
        lines=ifile.readlines()
        for line in lines:
            line=line.strip('\n')
            ID=line.split('\t',)[0]
            contig=line.split('\t',)[1]
            sstart=int(line.split('\t',)[8])
            send=int(line.split('\t',)[9])
            start=int(line.split('\t',)[6])
            end=int(line.split('\t',)[7])
            length=end-start
            if sstart < send:
                monomer=seq[contig][(sstart-1):send]
            else:
                monomer=seq[contig][(send-1):sstart]
            ofile.write(">%s.%s.%s-%s.%s\n%s\n" %(ID,contig,str(start),str(end),length,monomer))

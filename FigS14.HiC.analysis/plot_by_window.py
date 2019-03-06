#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This exript takes in .summary file from count_interaction_by_window.py, and output "summary_50kbwindow.txt" file for plotting in R.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

context_file="dmel_scaffold2_plus0310_context.txt"              
sizefile="dmel_scaffold2_plus0310.sizes"
path="Sexton_"
cens=["3R_5","Contig119","Contig79","Y_Contig26",'tig00057289']
outfile="summary.txt"

context={}
with open (context_file,'r') as cfile:    
        lines=cfile.readlines()
        for line in lines:
            line=line.strip('\n')
            if not line.startswith('#'):
                info=line.split('\t')              
                chro=info[0]
                het=info[1]
                region=info[2]
                name=info[3] 
                #info[4] is the range if it applied, some contigs are divided into two parts                
                if info[4]=='':
                    context[name]=(chro,het,region)
                # if the contig is divided into different parts
                else: 
                    if name not in context.keys():
                        context[name]={} 
                    context[name][info[4]]=(chro,het,region)
                    
with open (outfile,'w') as ofile:
    ofile.write("%s\t%s\t%s\t%s\t%s\t%s\n" %('#category','#chromosome','#contig_1','#contig_2','#window','#interactions per 100kb'))
    for cen in cens:

        belong=''
        if cen=='tig00057289' or cen=='Contig142':
            belong='centromere_2'
        elif cen=='3R_5' or cen=='tig00022795':
            belong='centromere_3'
        elif cen=='Contig119':
            belong='centromere_4'
        elif cen=='Contig79':
            belong='centromere_X'
        elif cen=='Y_Contig26':
            belong='centromere_Y'
            
        # define the context feature of the cen contig
        if cen!='Contig142':
            cen_chro=context[cen][0]
        else:
            cen_chro=context[cen]['0-60000'][0]
        #the path to contact files generated from HIC output 
        contact_file=path+cen+'.summary'
        with open(contact_file,'r') as ifile:
            lines=ifile.readlines()
            for line in lines:
                if not line.startswith('#'):
                    line=line.strip('\n')
                    info=line.split('\t')
                    contig=info[0].strip(' ')
                    window=info[1].strip(' ')
                    count=float(info[3].strip(' '))
                    # whole contig is one entry in context file
                    if type(context[contig])==tuple:
                        chro=context[contig][0]
                        het=context[contig][1]
                        region=context[contig][2]
                    # contig is divided into defferent regions, need to decide which region the window belongs to
                    elif type(context[contig])==dict:
                        window_s=int(window.split('-')[0])                       
                        window_e=int(window.split('-')[1])
                        length=window_e-window_s
                        for key in context[contig].keys():
                            key_s=int(key.split('-')[0])
                            key_e=int(key.split('-')[1])
                            # window totally overlap with key
                            if (window_s >=key_s and window_e <= key_e) or (window_s <= key_s and window_e >= key_e): 
                                chro=context[contig][key][0]
                                het=context[contig][key][1]
                                region=context[contig][key][2]
                                break
                            # window belongs to the contig with more than half of it overlaps
                            elif (key_s <= window_s <= key_e-length*0.5) or (key_s+length*0.5 <= window_e <= key_e):
                                chro=context[contig][key][0]
                                het=context[contig][key][1]
                                region=context[contig][key][2]
                                break
                        else:
                            print ("Error with the context file, can't locate"+contig+":"+window+".")
                    else:
                        print ("Error with the context file, can't locate"+contig+":"+window+".")
                    
                    #define category for each interaction
                    cat=''                    
                    
                    #don't count contigs without assigned chromosme
                    if chro=='':
                        pass  
                    
                    #don't count itself, or for Contig142, don't count the first 100 or 50kb                                                       
                    elif chro == cen_chro and region == 'centromere':
                        pass
                    
                    else:
                        # they are in the different chromosomes and it's in heterochromatin
                        if chro != cen_chro and het == 'het':                            
                            # and is also cen
                            if region=='centromere':
                                cat='centromere'
                            else:
                                cat='inter-heterochromatin' 
                        # they are in the same chromosome and in heterochromatin
                        elif chro == cen_chro and het == 'het':
                            # and is in proximal_centromeric region
                            if region == 'proximal_centromeric region':
                                cat='proximal_centromeric'
                            # and is in distal_centromeric region
                            else:
                                cat='distal_centromeric'
                        # inter- or intra- euchromatin regions
                        else:
                            cat='intra/inter-euchromatin'                        
                        
                        ofile.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(cat, belong, cen, contig, window,count))
    
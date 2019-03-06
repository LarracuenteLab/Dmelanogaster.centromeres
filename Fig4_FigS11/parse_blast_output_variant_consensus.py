#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script parse blast output from blasting repeat variants to consensus.

@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""


genome='dmel_scaffold2_plus0310.fasta'
gff='dmel_scaffold2_plus0310.repeat.gff3'
cens=["3R_5","Contig119","Contig79","Y_Contig26","tig00057289"]

file='G2'
consensus='G2_DM'

#change based on repeat full length
if consensus=='G2_DM':
    fullength='264-2957'
elif consensus=='Jockey3_Dsim':
    fullength='2148-4145'
elif consensus=='G_DM':
    fullength='1550-3782'
elif consensus=='G6_DM':
    fullength='322-1611'
elif consensus=='DOC2_DM':
    fullength='2029-4216'
elif consensus=='R1_DM':
    fullength='1728-4790'
elif consensus=='TART_A':
    fullength='6885-10124'
elif consensus=='HETA':
    fullength='1-6081'
elif consensus=='TAHRE':
    fullength='3386-6694' 
elif consensus=='PROTOP':
    fullength='1221-3785'
elif consensus=='BARI_DM':
    fullength='379-1398'        
elif consensus=='Jockey1_DSi':
    fullength='95-2629'
elif consensus=='DM1731_I':
    fullength='95-3812'
else:
    print ("unknown repeat")    


"""
#chcek both see if agree
#record variants' alignment coordinates along the consensus sequences
"""

blast=file+'.blast'

# one method: record alignment coordinates on query (check order and space between fragments see if combine)
summary={}
with open(blast,'r') as bfile:
    lines=bfile.readlines()
    for line in lines:
        info=line.split('\t',)
        variant=info[1]
        identity=int(float(info[2]))
        length=info[3]
        #aligment on consensus
        refstart=info[6]
        refend=info[7]
        #aligment on variant
        v_start=info[8]
        v_end=info[9]
        # mark the aligment order
        if int(v_start)<int(v_end):
            mark='f'
        else:
            mark='r'            
        if not variant in summary.keys():
            summary[variant]=(refstart+'.'+refend+'.'+length+'.'+str(identity))
            old_mark=mark
        else:            
            oldstart=summary[variant].split('.',)[0]
            oldend=summary[variant].split('.',)[1]
            #check if the same order 
            if mark == old_mark :
                # and with in 150bp
                if int(refend) > (int(oldstart)-150) and int(refstart) < int(oldstart): 
                    start=refstart
                    end=oldend
                elif int(refstart) < (150+int(oldend)) and int(refend) > int(oldend):
                    start=oldstart
                    end=refend                   
                else:
                    start=oldstart
                    end=oldend                    
                length=str(int(end)-int(start)+1)
                summary[variant]=(start+'.'+end+'.'+length+'.'+str(identity))

#method two: record alignment coordinates on query (combine all fragments)  
summary={}
with open(blast,'r') as bfile:
    lines=bfile.readlines()
    for line in lines:
        info=line.split('\t',)
        variant=info[1]
        identity=int(float(info[2]))
        length=info[3]
        refstart=info[6]
        refend=info[7]
        if not variant in summary.keys():
            summary[variant]=(refstart+'.'+refend+'.'+length+'.'+str(identity))
        else:
            old_identity=int(summary[variant].split('.',)[3])
            oldstart=summary[variant].split('.',)[0]
            oldend=summary[variant].split('.',)[1]
            if int(refstart) < int(oldstart):
                start=refstart
            else:
                start=oldstart
            if int(refend) > int(oldend):
                end=refend
            else:
                end=oldend
            length=str(int(end)-int(start)+1)
            
            if identity < old_identity:
                new_identity=identity
            else:
                new_identity=old_identity
            summary[variant]=(start+'.'+end+'.'+length+'.'+str(new_identity))
    
"""
output summary file from blast output
"""
#summary blast output (note that summary dictionary records alignment coordinates on query(consensus))

sumfile=file+'_blast.summary'

test=[]
with open(sumfile,'w') as sfile:
    sfile.write("#%s\t#%s\t#%s\t#%s\n" %("variant","contig-start-end-length","consensus-start-end-length","identity%"))
    for variant in summary.keys():
        info=summary[variant].split('.',)
        # alignment on consensus
        consensus_start=info[0]
        consensus_end=info[1]
        consensus_length=info[2]
        identity=info[3]
        #only output sequences with > 80% identical to consensus
        if int(identity) < 80:
            print (variant+': low identity '+identity)
            test.append(variant)
        else:
            contig=variant.split('.',)[1]
            genome_start=variant.split('.',)[2].split('-',)[0]
            genome_end=variant.split('.',)[2].split('-',)[1]
            genome_length=variant.split('.',)[3]
            sfile.write("%s\t%s-%s-%s-%s\t%s-%s-%s-%s\t%s\n" %(variant,contig,genome_start,genome_end,genome_length,consensus,consensus_start,consensus_end,consensus_length,identity))


# if miss annotation, output sequences
sequence=file+'_variance.fasta'
miss=file+'_missannotation.fasta'
with open(miss,'w') as ofile:           
    with open(sequence,'r') as ifile:
        lines=ifile.readlines()
        for line in lines:
            if line.startswith('>'):
                line=line.strip('\n>')
                if line in test:
                    i=lines.index('>'+line+'\n')
                    seq=lines[i+1].strip('\n')
                    ofile.write(">%s\n%s\n" %(line,seq))

 
"""
output sumamry file for plotting
"""         
# output location information for plotting (note sumfile contain info of "#variant,contig-start-end-length,consensus-start-end-length,identity%")

plotfile=file+'_plot.summary'

full_count=0
total_count=0
with open(plotfile,'w') as ofile:
    ofile.write("#%s\t#%s\t#%s\t#%s\n" %("name","contig","location","length"))
    with open(sumfile,'r') as sfile:
        lines=sfile.readlines()
        for line in lines:
            if not line.startswith('#'):
                line=line.strip('\n')
                info=line.split('\t')
                #TE name
                name=info[0].split('.',)[0]
                #location info: contig, start position
                contig=info[1].split('-',)[0]                
                location=info[1].split('-',)[1]
                # define if it's full length
                consensus_start=info[2].split('-',)[1]
                consensus_end=info[2].split('-',)[2]
                if int(consensus_start) <= start and int(consensus_end) >= end:
                    full='full'
                    full_count+=1
                elif int(consensus_length) > 1000:
                    full='big'
                else:
                    full='small'
                total_count+=1
                ofile.write("%s\t%s\t%s\t%s\n" %(name,contig,location,full))

print ("total count: "+str(total_count)+"\n")
print ("full length: "+str(full_count)+"\n")


"""
category TE locations based on context file
output statistics summary(count of total and full-length elements)
"""

context_file="dmel_scaffold2_plus0310_context.txt"

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
                leng=info[5]
                #info[4] is the range if it applied, some contigs are divided into two parts                
                if info[4]=='':
                    context[name]=(chro,het,region,leng)
                # if the contig is divided into different parts
                else:
                    leng=int(info[4].split('-',)[1])-int(info[4].split('-',)[0])
                    if name not in context.keys():
                        context[name]={}                    
                    context[name][info[4]]=(chro,het,region,leng)

# get size of het/eu chromatin regions in genome
size={}
size['eu']=0
size['het']=0
size['centromere']=0
size['telomere']=0
for key in context.keys():
    #contig is divided into regions
    if type(context[key])==dict:
        for region in context[key].keys():
            #euchromatin
            if context[key][region][1]=='eu':
                size['eu']+=int(context[key][region][3])
            #heteorchromatin
            elif context[key][region][1]=='het':
                size['het']+=int(context[key][region][3])    
                if context[key][region][2]=='centromere':
                    size['centromere']+=int(context[key][region][3])
                elif context[key][region][2]=='telomere':                        
                    size['telomere']+=int(context[key][region][3])
            else:
                print ("unknown contig: "+key)
    #contig is not divided into regions
    else:
        if context[key][1]=='eu':
            size['eu']+=int(context[key][3])
        #heteorchromatin
        elif context[key][1]=='het':
            size['het']+=int(context[key][3])    
            if context[key][2]=='centromere':                        
                size['centromere']+=int(context[key][3])
            elif context[key][2]=='telomere':                        
                size['telomere']+=int(context[key][3])
        else:
            print ("unknown contig: "+key)
            
size_file="contig_size.txt"
with open(size_file,'w') as sfile:
    for key in size.keys():
        sfile.write("%s\t%s\n" %(key,str(size[key])))
    

plotfile=file+'_plot.summary'
staticsfile=file+'.statistics'

total_count=0
count={}
count['eu']={}
count['pericentromeric-heterochromatin']={}
count['centromere']={}
count['telomere']={}

for key in count.keys():
    count[key]['full']=0
    count[key]['truncated']=0    

with open(staticsfile,'w') as s_file:
    with open(plotfile,'r') as p_file:
        lines=p_file.readlines()
        for line in lines:
            if not line.startswith('#'):
            #if not line.startswith('#') and (not line.split('\t',)[1].startswith('Y')):
                cat=''
                line=line.strip('\n')
                info=line.split('\t')
                name=info[0]
                contig=info[1]
                pos=info[2]
                full=info[3]
                
                #total count
                total_count+=1
    
                # define category (eu,het,centromere,telomere)
                for key in context.keys():
                    if key==contig:
                        #contig is divided into regions
                        if type(context[key])==dict:
                            for region in context[key].keys():
                                start=region.split('-',)[0]
                                end=region.split('-',)[1]
                                if int(start) <= int(pos) <= int(end):
                                    category=context[key][region]
                        else:
                            category=context[key]
                        break
                chro=category[0]
                if category[1]=='eu':
                    cat='eu'
                elif category[1]=='het':
                    if category[2]=='centromere':                        
                        cat='centromere'
                    elif category[2]=='telomere':                        
                        cat='telomere'
                    else:
                        cat='pericentromeric-heterochromatin'
                else:
                    print ('unknow eu/het')
                
                # count each cat
                if full=='full':
                    count[cat]['full']+=1
                else:
                    count[cat]['truncated']+=1
                
                s_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(name,chro,cat,contig,pos,full))
    
    s_file.write("\n#%s: %s\n" %(consensus,total_count))
    for key in count.keys():
        for item in count[key].keys():
            s_file.write("#%-30s\t%-10s\t%s\n" %(key,item,str(count[key][item])))





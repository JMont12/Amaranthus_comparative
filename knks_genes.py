#!/usr/bin/python

#this script will take in three tsv files output by Synmap, calculate kn/ks values and output useful information about how many genes have kn/ks values >1 and are common amoungst species
#this script will output the number of genes to make a ven diagram to std out. It will also write the gene names of all genes in common amongst the three species to the provided out_file path.
#in_files should be formatted in the way that CoGe Synmap outputs in the "Results with synonymous/non-synonymous rate values" link
#usage /path/to/knks_genes.py in_SP.txt in_WH.txt in_PA.txt /path/hypo_gff /path/to/out.txt

#import some important packages
from sys import argv, version_info
from os.path import realpath, splitext
import os

#this line parses the command line and stores the file paths that you give as variables that can be called later
in_1, in_2, in_3, gff_in, out_file = realpath(argv[1]), realpath(argv[2]), realpath(argv[3]), realpath(argv[4]), realpath(argv[5])

fh_1=open(in_1, 'r')
fh_2=open(in_2, 'r')
fh_3=open(in_3, 'r')
fh_gff=open(gff_in, 'r')
out_fh=open(out_file, 'w+')

#this dictionary will store all genes with a knks >1 in any of the three species
#to where each gene name is the key, the value will be an array with [gene name, SP, WH, PA presence/absence (1/0)]
genes={}

gene_names={}

SP_count=0
WH_count=0
PA_count=0
SP_gene_count=0
WH_gene_count=0
PA_gene_count=0
SP_undefined_count=0
WH_undefined_count=0
PA_undefined_count=0
SP_same_count=0
WH_same_count=0
PA_same_count=0
fake_lines=0
line_count=1

for line in fh_gff:
        if '\tgene\t' in line:
                line=line.strip('\n')
                parts=line.split('\t')
                subparts=parts[8].split('=')
                name_parts=subparts[1].split(';')
                if len(subparts) >7:
	                gene_names[str(name_parts[0])+'-RA']=subparts[7]
			if 'RB' in name_parts[0]:
				gene_names[str(name_parts[0])+'-RB']=subparts[7]
		else:
			gene_names[str(name_parts[0])+'-RA']='Protein of unknown function'
                        if 'RB' in name_parts[0]:
                                gene_names[str(name_parts[0])+'-RB']=subparts[7]
fh_gff.close

#hardcode a value for a gene identified in the palmer amaranth analysis

gene_names['AH007694-RB']='Protein of unknown function'


for line in fh_1:
	if line_count>3:
		parts=line.split('\t')
		subparts=parts[3].split('|')
		try:
			float(parts[0])
			SP_gene_count+=1
                        ks=float(parts[0])
                        kn=float(parts[1])
                        if ks>0.0:
				if kn/ks >1:
					if subparts[6] not in genes:
						genes[subparts[6]]=[subparts[6],1,0,0]
					else:
						genes[subparts[6]][1]+=1
                               		SP_count+=1
               		else:
                        	if float(parts[1])>0.0000:
					if subparts[6] not in genes:
						genes[subparts[6]]=[subparts[6],1,0,0]
					else:
						genes[subparts[6]][1]+=1

                                	SP_undefined_count+=1
				else:
					SP_same_count+=1
		except ValueError:
			fake_lines+=1

	line_count+=1
fh_1.close

print "done with hybridus"

line_count=1

for line in fh_2:
        if line_count>3:
                parts=line.split('\t')
                subparts=parts[3].split('|')
                try:
                        float(parts[0])
                        WH_gene_count+=1
                        ks=float(parts[0])
                        kn=float(parts[1])
                        if ks>0.0:
                                if kn/ks >1:
	                                if subparts[6] not in genes:
        	                                genes[subparts[6]]=[subparts[6],0,1,0]
						WH_count+=1
                                	else:
						genes[subparts[6]][2]+=1
						WH_count+=1
                        else:
                                if float(parts[1])>0.0000:
                                        if subparts[6] not in genes:
                                                genes[subparts[6]]=[subparts[6],0,1,0]
                                                WH_undefined_count+=1
                                        else:
                                                genes[subparts[6]][2]+=1
						WH_undefined_count+=1
				else:	
					WH_same_count+=1
                except ValueError:
                        fake_lines+=1

        line_count+=1
fh_2.close

print "done with tuberculatus"
line_count=1

for line in fh_3:
        if line_count>3:
                parts=line.split('\t')
                subparts=parts[3].split('|')
                try:
                        float(parts[0])
                        PA_gene_count+=1
                        ks=float(parts[0])
                        kn=float(parts[1])
                        if ks>0.0:
                                if kn/ks >1:
                                        if subparts[6] not in genes:
                                                genes[subparts[6]]=[subparts[6],0,0,1]
                                                PA_count+=1
                                        else:
                                                genes[subparts[6]][3]+=1
                                                PA_count+=1
                        else:
                                if float(parts[1])>0.0000:
                                        if subparts[6] not in genes:
                                                genes[subparts[6]]=[subparts[6],0,0,1]
                                                PA_undefined_count+=1
                                        else:
                                                genes[subparts[6]][3]+=1
						PA_undefined_count+=1
                                else:
                                        PA_same_count+=1
                except ValueError:
                        fake_lines+=1
        line_count+=1
fh_3.close

SP_only=0
WH_only=0
PA_only=0
SP_WH=0
SP_PA=0
PA_WH=0
all=0

out_fh.write("hypochondriacus_gene_name"+"\t"+"hybridus_presence"+"\t"+"tuberculatus_presence"+"\t"+"palmeri_presence"+"\n")
for gene in genes:
	out_fh.write(str(genes[gene][0])+"\t"+str(genes[gene][1])+"\t"+str(genes[gene][2])+"\t"+str(genes[gene][3])+"\t"+str(gene_names[gene])+"\n")
	if genes[gene][1]>0 and genes[gene][2]==0 and  genes[gene][3]==0:
		SP_only+=1
	if genes[gene][1]==0 and genes[gene][2]>0 and  genes[gene][3]==0:
		WH_only+=1
	if genes[gene][1]==0 and genes[gene][2]==0 and  genes[gene][3]>0:
		PA_only+=1
	if genes[gene][1]>0 and genes[gene][2]>0 and  genes[gene][3]==0:
		SP_WH+=1
	if genes[gene][1]>0 and genes[gene][2]==0 and  genes[gene][3]>0:
		SP_PA+=1
	if genes[gene][1]==0 and genes[gene][2]>0 and  genes[gene][3]>0:
		PA_WH+=1
	if genes[gene][1]>0 and genes[gene][2]>0 and  genes[gene][3]>0:
		all+=1

print "there were "+str(SP_count)+" smooth pigweed genes under positive selection out of "+str(SP_gene_count)+" gene pairs, with "+str(SP_undefined_count)+" more having no synonymous mutations but some non-synonymous. "+str(SP_same_count)+" genes were identical to hypochondriacus"
print "there were "+str(WH_count)+" waterhemp genes under positive selection out of "+str(WH_gene_count)+" gene pairs, with "+str(WH_undefined_count)+" more having no synonymous mutations but some non-synonymous. "+str(WH_same_count)+" genes were identical to hypochondriacus"	
print "there were "+str(PA_count)+" palmer amaranth genes under positive selection out of "+str(PA_gene_count)+" gene pairs, with "+str(PA_undefined_count)+" more having no synonymous mutations but some non-synonymous. "+str(PA_same_count)+" genes were identical to hypochondriacus"
print ""
print "ven diagram information:"
print ""
print "There were "+str(SP_only)+" genes specific to smooth pigweed"
print "There were "+str(WH_only)+" genes specific to waterhemp"
print "There were "+str(PA_only)+" genes specific to palmer amaranth"
print "There were "+str(SP_WH)+" genes specific to smooth pigweed and waterhemp"
print "There were "+str(SP_PA)+" genes specific to smooth pigweed and palmer amaranth"
print "There were "+str(PA_WH)+" genes specific to waterhemp and palmer amaranth"
print "There were "+str(all)+" genes shared between all species"

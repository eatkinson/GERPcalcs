### calculate the average GERP score in UCSC canonical transcript genes. Do for all sites in the HGDP dataset annotated with their RS scores. Has the whole genome so nicer than chr by chr which is how 1000G is organized. Probably fine to use these, will have fewer variants than 1000G but good enough.

#load in all the modules I am going to need
import numpy as np
import vcf
import string
import re

outFile = open('hg19_avgGERPintergenic_HGDP.txt','w')

#load in the HGDP genomes Gerp, or really RS-annotated 
vcf = vcf.Reader(open('all_combined.snps.autos_and_PAR.vqsr99.BEAGLE.snpEff.RefAnc.GERP.callableMask.vcf.gz','r'))
#and load in the bed file of windows - here, the start and stop positions of the canonical transcripts 
winfile = open("hg19_intergenic.bed", "r")
wins = np.genfromtxt(winfile, dtype=None)  #import the text of the windows file as a numpy array


for line in range(0, len(wins)):
#look up SNPs in the VCF file in terms of the windows. The first col is the transript ID, then the chrom and then start and stop pos for the canonical transcript
    chrom = wins[line][0]
    start = wins[line][1]
    stop = wins[line][2]
    
    count = 0  #count is the number of SNPs in this window
    gerpavg = 0.0
    gerpsum = 0.0  #ensure this will be output in float format
    
    for record in vcf.fetch(chrom,start,stop):  #get all the SNPs in the vcf file that are in this window
        count = count + 1
    	#gerpsum = 0.0 
        gerpsum = gerpsum + float(record.INFO['RS'][0])
    if count > 0 :
    	gerpavg = gerpsum/count
    else : gerpavg = "NA"
    output = "%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(start), str(stop), str(gerpavg), str(count))
    #outputs the chromosome, the start and stop sites of the particular window, the average GERP score for SNPs in that window, 
        
    outFile.write(output)
outFile.close()

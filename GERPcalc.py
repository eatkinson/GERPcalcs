#python script to calculate GERP score averages for 600bp windows

import numpy as np
import vcf
import string
import re

vcf = vcf.Reader(open('/vault/henn/genomes/hgdp/variants/all_combined.snps.autos_and_PAR.vqsr99.BEAGLE.snpEff.RefAnc.GERP.callableMask.vcf.gz','r'))
outFile = open('GERPtest.txt','w')
winfile = open("/vault/henn/people/elizabeth/FOXP2/Phasing_Network/chimps/HGDPgenomes_GERP.INFO.Windows.txt", "r")

wins = np.genfromtxt(winfile, dtype=None)  #import the text of the windows file as a numpy array

for line in range(0, len(wins)): #look up SNPs in the VCF file in terms of the windows, since these are overlapping sliding windows.
    
    start = wins[line][1]
    stop = wins[line][2]
    chrom = wins[line][0]
    
    count = 0  #count is the number of SNPs in this window
    HighRS = 0   #HighRS is the number of SNPs in this window that have a high GERP score, i.e. >= 3
    
    for record in vcf.fetch(chrom,start,stop): #get all the SNPs in the vcf file that are in this window
        count = count + 1    
        gerpsum = 0.0  #ensure this will be output in float format
        gerpsum = gerpsum + float(record.INFO['RS'][0])
        if record.INFO['RS'][0] >= 3:
            HighRS = HighRS + 1
    gerpavg = gerpsum/count
    output = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(start), str(stop), str(gerpavg), str(count), str(HighRS))
    #outputs the chromosome, the start and stop sites of the particular window, the average GERP score for SNPs in that window, how many SNPs are in that window, and how many SNPs have GERP scores over 3 in that window.
    
    outFile.write(output)
outFile.close()
  
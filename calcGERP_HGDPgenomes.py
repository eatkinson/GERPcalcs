
# coding: utf-8

# In[9]:


import numpy as np
import string
import re


# In[18]:


GERPinfo = open('HGDPgenomes_noncoding_maf5_het_RS_2.txt','r')
    #GERPinfo has 5 columns: the chromosome, position, heterozygosity of the SNP,
    #alleles for REF and ALT, and RS score of the SNPs individually
GERPs = np.genfromtxt(GERPinfo, dtype=None)   #turn the GERP files into numpy arrays


# And then can use the window file to calculate the average GERP, avg and sd heterozygosity, and count of SNPs in each window. These SNPs are already filtered to only be in non-coding regions.
# load in the bed file of windows 
winfile = open("HGDPgenomes_2251bpWins.txt", "r")
wins = np.genfromtxt(winfile, dtype=None)  #import the text of the windows file as a numpy array


# In[ ]:


outFile = open('avgGERP_HGDPnoncoding_2251wins_notchr1.txt','w')

for line in range(0, len(wins)):

    start = wins[line][1]
    stop = wins[line][2]
    chrom = wins[line][0]

    count = 0  #count is the number of SNPs in this window
    HighRS = 0  #HighRS is the number of SNPs in this window that have a high GERP score, i.e. >= 3

    for item in range(0, len(GERPs)):
        pos = GERPs[item][1]
        RS = GERPs[item][5]
        GERPchr = GERPs[item][0]
        het =  GERPs[item][2]
        
        if chrom != "chr1" :
            if chrom == GERPchr:
                if pos in range(start, stop) :
                    count = count + 1
                    hetsum = 0.0
                    gerpsum = 0.0  #ensure this will be output in float format
                    gerpsum = gerpsum + RS
                    hetsum = hetsum + het

                    if RS >= 3:
                        HighRS = HighRS + 1

    if count > 0 :
        gerpavg = gerpsum/float(count)
        hetavg = hetsum/float(count)
        
    else:
        gerpavg = "NA"
        hetavg = "NA"

    output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(start), str(stop), str(count), str(HighRS),str(gerpavg),str(hetavg))
    #outputs the chromosome, the start and stop sites of the particular window, the average GERP score for SNPs in that window, 
    #how many SNPs are in that window, and how many SNPs have GERP scores over 3 in that window.
    #currently tailored to only output windows where there is at least 1 SNP

    outFile.write(output)
outFile.close()



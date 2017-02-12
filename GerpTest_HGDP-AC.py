
# coding: utf-8

# In[37]:

get_ipython().system(u'pwd')


# In[44]:

import numpy as np
import vcf
import string
import re


# In[45]:

#HGDP genomes VCF file
vcf = vcf.Reader(open('HGDPgenomes_GeneMask.recode.vcf.gz','r'))
#this is the HGDP genomes file that I have filtered to only have the genic areas in it.


# In[46]:

outFile= open('HGDPgenomes_GERPtest2.txt', 'w')


# In[47]:

#HGDP genomes file:
winfile = open("/vault/henn/people/elizabeth/FOXP2/SelectionTests/GERP/HGDPgenomes_filteredWindows", "r")
#winfile should have columns of the chromosome, start position, and stop position of all of the windows


# In[ ]:

wins = np.genfromtxt(winfile, dtype=None)  #import the text of the windows file as a numpy array


# In[ ]:

for line in range(0, len(wins)):
#look up SNPs in the VCF file in terms of the windows, since these are overlapping sliding windows.
    
    start = wins[line][1]
    stop = wins[line][2]
    chrom = wins[line][0]
    
    count = 0  #count is the number of SNPs in this window
    HighRS = 0  #HighRS is the number of SNPs in this window that have a high GERP score, i.e. >= 3
    
    for record in vcf.fetch(chrom,start,stop):  #get all the SNPs in the vcf file that are in this window
        RS = float(record.INFO['RS'][0])   #RS is the annotation of the SNP's RS score
        AC = float(record.INFO['AC'][0])   #AC is the allele count for the alternate allele at the site. 
        
        if AC > 1:   #only count SNPs that are not singletons. Those could be one-off sequencing errors
            count = count + 1   #tally up the SNPs in the window
            gerpsum = 0.0  #ensure this will be output in float format
            gerpsum = gerpsum + RS
            
            if RS >= 3:
                HighRS = HighRS + 1   #tally up the high RS scoring SNPs
    
    if count > 0 :
        gerpavg = gerpsum/float(count)
   
        output = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(start), str(stop), str(gerpavg), str(count), str(HighRS))
        #outputs the chromosome, the start and stop sites of the particular window, the average GERP score for SNPs in that window, 
        #how many SNPs are in that window, and how many SNPs have GERP scores over 3 in that window.
        
    outFile.write(output)
outFile.close()
  


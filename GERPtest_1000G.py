import numpy as np

winfile = open('/vault/henn/people/elizabeth/FOXP2/SelectionTests/GERP/WinFile.chr1', 'r') #%d is the dictionary character
wins = np.genfromtxt(winfile, dtype=None)  #import the text of the windows file as a numpy array
gerpFile1 = open('/vault/henn/people/elizabeth/FOXP2/SelectionTests/GERP/1000genomes_chr1.GERPS.txt', 'r')
gerps = np.genfromtxt(gerpFile1, dtype=None)
outFile = open('GERP_test1', 'w')

for line in range(0, len(wins)):

    start = wins[line][1]
    stop = wins[line][2]
    chrom = wins[line][0]

    count = 0  #count is the number of SNPs in this window
    HighRS = 0  #HighRS is the number of SNPs in this window that have a high GERP score, i.e. >= 3

    for item in range(0, len(gerps)):
        pos = gerps[item][1]
        RS = gerps[item][2]
        if pos in range(start, stop) :
            count = count + 1
            gerpsum = 0.0  #ensure this will be output in float format
            gerpsum = gerpsum + RS

            if RS >= 3:
                HighRS = HighRS + 1

    if count > 0 :
        gerpavg = gerpsum/float(count)
    else:
        gerpavg = "NA"

    output = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chrom), str(start), str(stop), str(gerpavg), str(count), str(HighRS))
    #outputs the chromosome, the start and stop sites of the particular window, the average GERP score for SNPs in that window, 
    #how many SNPs are in that window, and how many SNPs have GERP scores over 3 in that window.
    #currently tailored to only output windows where there is at least 1 SNP

    outFile.write(output)
outFile.close()


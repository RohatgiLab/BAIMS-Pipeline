"""
Developed by Bhaven Patel, Rajat Rohatgi Lab
This code is free software; you can redistribute it and/or modify it.
Please cite when using or altering this code.
"""
#This script generates general statistics about alignments and mismatches for alignments from a .sam file.

import sys, getopt

# chromosomeCount function determines which chromosome the read aligns to
def chromosomeCount(tokens, chromList, extraChromList):
    chromosome = tokens[2];
    if tokens[2] == "chrM" or len(tokens[2])>5:  #for mitochondrial DNA or unassembled region
        for x in range(0,len(extraChromList),1):
            if extraChromList[x] == tokens[2]: #if region has been seen before in alignments  
                chromList[23+1+x]+=1;
                return;
        extraChromList.append(tokens[2]); #if region has not been seen before in alignments                 
        chromList.append(1);
    else:
        if chromosome[3:] == "X":
            chromList[22]+=1;
        elif chromosome[3:] == "Y":
            chromList[23]+=1;
        else:
            num = int(chromosome[3:]);
            chromList[num-1]+=1;

# writeSummaryStats function writes stats into output file
def writeSummaryStats(output, alignedList, multiAlignments, chromList, extraChromList, uniqueAlignments):
    output.write("Number of reads:\t\t\t" +str(alignedList[0]+alignedList[1])+"\n");
    output.write("Number of aligned reads:\t\t\t" +str(alignedList[1])+"\n");
    output.write("Number of reads that did not align:\t\t\t" +str(alignedList[0])+"\n");
    output.write("Number of reads aligned >1 times:\t\t\t" +str(multiAlignments)+"\n");
    output.write("Number of reads aligned uniquely:\t\t\t" +str(uniqueAlignments[0]+uniqueAlignments[1]+uniqueAlignments[2]+uniqueAlignments[3])+"\n");
    output.write("Number of unique alignments with 0 mismatches:\t\t\t" +str(uniqueAlignments[0])+"\n");
    output.write("Number of unique alignments with 1 mismatches:\t\t\t" +str(uniqueAlignments[1])+"\n");
    output.write("Number of unique alignments with 2 mismatches:\t\t\t" +str(uniqueAlignments[2])+"\n");
    output.write("Number of unique alignments with 3 mismatches:\t\t\t" +str(uniqueAlignments[3])+"\n\n");
    output.write("Number of unique alignments per chromosome:\n");
    output.write("Number of reads aligned to chrom1:\t\t\t" +str(chromList[0])+"\n");
    output.write("Number of reads aligned to chrom2:\t\t\t" +str(chromList[1])+"\n");
    output.write("Number of reads aligned to chrom3:\t\t\t" +str(chromList[2])+"\n");
    output.write("Number of reads aligned to chrom4:\t\t\t" +str(chromList[3])+"\n");
    output.write("Number of reads aligned to chrom5:\t\t\t" +str(chromList[4])+"\n");
    output.write("Number of reads aligned to chrom6:\t\t\t" +str(chromList[5])+"\n");
    output.write("Number of reads aligned to chrom7:\t\t\t" +str(chromList[6])+"\n");
    output.write("Number of reads aligned to chrom8:\t\t\t" +str(chromList[7])+"\n");
    output.write("Number of reads aligned to chrom9:\t\t\t" +str(chromList[8])+"\n");
    output.write("Number of reads aligned to chrom10:\t\t\t" +str(chromList[9])+"\n");
    output.write("Number of reads aligned to chrom11:\t\t\t" +str(chromList[10])+"\n");
    output.write("Number of reads aligned to chrom12:\t\t\t" +str(chromList[11])+"\n");
    output.write("Number of reads aligned to chrom13:\t\t\t" +str(chromList[12])+"\n");
    output.write("Number of reads aligned to chrom14:\t\t\t" +str(chromList[13])+"\n");
    output.write("Number of reads aligned to chrom15:\t\t\t" +str(chromList[14])+"\n");
    output.write("Number of reads aligned to chrom16:\t\t\t" +str(chromList[15])+"\n");
    output.write("Number of reads aligned to chrom17:\t\t\t" +str(chromList[16])+"\n");
    output.write("Number of reads aligned to chrom18:\t\t\t" +str(chromList[17])+"\n");
    output.write("Number of reads aligned to chrom19:\t\t\t" +str(chromList[18])+"\n");
    output.write("Number of reads aligned to chrom20:\t\t\t" +str(chromList[19])+"\n");
    output.write("Number of reads aligned to chrom21:\t\t\t" +str(chromList[20])+"\n");
    output.write("Number of reads aligned to chrom22:\t\t\t" +str(chromList[21])+"\n");
    output.write("Number of reads aligned to chromX:\t\t\t" +str(chromList[22])+"\n");
    output.write("Number of reads aligned to chromY:\t\t\t" +str(chromList[23])+"\n");
    if len(extraChromList) > 0: #print out anything in extraChromList                            
        increment = 0;
        for word in extraChromList:
            output.write("Number of reads aligned to " +word+":\t\t\t" +str(chromList[23+increment+1])+"\n");
            increment+=1;

               
# main method
def main(argv):
    #get arguments
    inputfile = '';
    outputfile= '';
    try:
        opts, args = getopt.getopt(argv,"hi:o:", ["ifile="]);
    except getopt.GetoptError:
        print 'Usage: analyzeSAM.py -i <inputfile> -o <outputfile>';
        sys.exit(2);
    for opt, value in opts:
        if opt == '-h':
            print 'analyzeSAM.py -i <inputfile> -o <outputfile>';
            sys.exit();
        elif opt in ("-i", "--ifile"):
            inputfile = value;
        elif opt in ("-o"):
            outputfile = value;
    print 'Input file is ', inputfile;
    #open inputfile
    input = open(inputfile,"r");
    alignedList = [0,0]; #holds number of reads that didn't-align/align; alignedList[0] = # reads that didn't align; alignedList[1] = # reads that aligned
    uniqueAlignments = [0]*4; #holds the number of uniquely aligned reads with 0,1,2, or 3 mismatches
    chromList = [0]*24; #array to hold the number of uniquely aligned reads per chromosome
    extraChromList = []; #array to hold the names of chromosome segments not fully integrated into the genome
    multiAlignments = 0;
    previous = ['string']; #variable to check if the current alignment comes from a read that aligns to multiple locations
    previous2 = ['string']; #variable to track if current alignment comes from a read that has already been counted as aligning to multiple locations
    line = input.readline(); #get first line of SAM file
    while line != "":
        tokens = line.split();
        # check if this alignment comes from a sequencing read that aligns to multiple locations in the genome since Bowtie should give the 5 best alignments; if it does, then move to next alignment in file
        if tokens[0] == previous[0]:
            previous = tokens;
            line = input.readline();
            continue;
        #alignment comes from a previously unseen sequencing read, so need to determine if it aligned uniquely    
        nextLine = None;
        # if-else counts number of reads that aligned
        if tokens[2] == "*":
            alignedList[0]+=1; #read didnt align
        else: #read did align, either uniquely or multiple times
            alignedList[1]+=1;
            #need to check if alignment is unique 
            nextLine = input.readline(); #get next line to see if it is an alignment from the same sequencing read
            if nextLine != "": 
                nextTokens = nextLine.split();
                if tokens[0] == nextTokens[0]: #check if current alignment and next alignment come from the same sequencing read
                    multiAlignments+=1; #increment counter for reads that align to more than one location
                else: #this alignment is unique so determine which chromosome it comes from
                    chromosomeCount(tokens, chromList,extraChromList);
                    #if-elif counts number of reads that aligned uniquely with 0,1,2, or 3 mismatches
                    if "NM:i:0" in tokens:
                        uniqueAlignments[0]+=1;
                    elif "NM:i:1" in tokens:
                        uniqueAlignments[1]+=1;
                    elif "NM:i:2" in tokens:
                        uniqueAlignments[2]+=1;
                    elif "NM:i:3" in tokens:
                        uniqueAlignments[3]+=1;
        #update line and previous check
        previous = tokens;
        if nextLine != None: #check if the nextLine was already retrieved
            line = nextLine;
        else:
            line = input.readline();

    input.close();

    # write statistics to output file
    output = open(outputfile, "w");
    writeSummaryStats(output,alignedList,multiAlignments,chromList,extraChromList,uniqueAlignments);
    output.close;

    print "Number of reads not aligned, aligned:",alignedList;
    print "Done with initial analysis of aligned sequencing reads";


#main is envoked    
main(sys.argv[1:]);

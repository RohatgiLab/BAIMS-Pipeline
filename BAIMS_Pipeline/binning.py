"""
Developed by Bhaven Patel, Rajat Rohatgi Lab
This code is free software; you can redistribute it and/or modify it.
Please cite when using or altering this code.
"""

import sys, getopt, re

# returns the proper bin index depending on the chromosome string passed in         
def getChromBinNum(chromosome, extraChromList, bins):
    chromNum = -1;
    if chromosome == "X": #to account for X and Y chromosomes              
        chromNum = 22;
    elif chromosome == "Y":
        chromNum = 23;
    elif chromosome == "M" or len(chromosome) > 2:
        for x in range(0,len(extraChromList),1):
            if extraChromList[x] == chromosome:
                return 23+x+1;
        extraChromList.append(chromosome);
        bins.append([]);
        return len(bins)-1;
    else:
        chromNum = int(chromosome)-1;
    return chromNum;       
     
# tallies up number of reads that align to each bin for each chromosome and keeps track of the number of reads with 0,1,2,3 mismatches per bin per chromosome
def addToBins(chromNum, startPos, bins, binSize, mismatches, multiAlign,sense):
    index = startPos/binSize; #calculates the index of the gene block the read should be put in
    if index >= len(bins[chromNum]): #adds new block to chromosome block list, if needed 
        for x in range(len(bins[chromNum])-1, index):
            mismatchesBlock = [];
            for  x in range(4): #creates list of unique vs non-unique alignment dependent on sense vs anti-sense alignments
                uniqueList = [];
                for x in range(2): #creates list of sense vs anti-sense insertions in the uniqueList
                    senseList = [0]*2;
                    uniqueList.append(senseList);
                mismatchesBlock.append(uniqueList);
            bins[chromNum].append(mismatchesBlock); 
    block = (bins[chromNum])[index];
    if multiAlign:
        if sense:
            block[mismatches][1][0]+=1; #increment non-unique sense count
        else:
            block[mismatches][1][1]+=1; #increment non-unique anti-sense count
    else:
        if sense:
            block[mismatches][0][0]+=1; #increment unique sense count
        else:
            block[mismatches][0][1]+=1; #increment unique anti-sense count


#prints bins when no offset specified
def printBinning(output, bins, binSize,extraChromList,annotationsDict):
    totalAlign = 0;
    unique = 0;
    for chromosome in range(0,len(bins),1):
        chromosomeStr = "";
        if chromosome == 22:
            chromosomeStr = "X";
        elif chromosome == 23:
            chromosomeStr = "Y";
        elif chromosome > 23:
            chromosomeStr = extraChromList[chromosome - 24];
        else:
            chromosomeStr = str(chromosome+1);
        chromosomeBlock = bins[chromosome];
        blockLen = len(chromosomeBlock);
        for cBlockInd in range(0, blockLen):
            output.write("Chromosome "+chromosomeStr+", "+str(cBlockInd*binSize)+"-"+str((cBlockInd+1)*binSize-1)+"bps:\n");
            output.write(("").ljust(15)+("Total alignments").ljust(20)+("+").ljust(10)+("-").ljust(10)+("Unique").ljust(15)+("+").ljust(10)+("-").ljust(10)+("Non-Unique").ljust(15)+("+").ljust(10)+("-").ljust(10)+"\n");
            mismatchesBlock = chromosomeBlock[cBlockInd];
            totalforBlock = 0;
            totalSense = 0;
            totalAnti = 0;
            totalUnique = 0;
            totalUnSense = 0;
            totalUnAnti = 0;
            totalNonUnq = 0;
            totalNonUnqSense = 0;
            totalNonUnqAnti = 0;
            for x in range(4):
                output.write((str(x)+" mismatches").ljust(15));
                output.write(str(mismatchesBlock[x][0][0]+mismatchesBlock[x][0][1]+mismatchesBlock[x][1][0]+mismatchesBlock[x][1][1]).ljust(20)); #total alignments                           
                output.write(str(mismatchesBlock[x][0][0]+mismatchesBlock[x][1][0]).ljust(10)); #total sense alignments                                                                       
                output.write(str(mismatchesBlock[x][0][1]+mismatchesBlock[x][1][1]).ljust(10)); #total antisense alignments                                                                   
                output.write(str(mismatchesBlock[x][0][0]+mismatchesBlock[x][0][1]).ljust(15)); #total unique alignments                                                                      
                output.write(str(mismatchesBlock[x][0][0]).ljust(10)); #total unique sense alignments                                                                                         
                output.write(str(mismatchesBlock[x][0][1]).ljust(10)); #total unique antisense alignments                                                                                     
                output.write(str(mismatchesBlock[x][1][0]+mismatchesBlock[x][1][1]).ljust(15)); #total non-unique alignments                                                                  
                output.write(str(mismatchesBlock[x][1][0]).ljust(10)); #total non-unique sense alignments                                                                                     
                output.write(str(mismatchesBlock[x][1][1]).ljust(10)+"\n"); #total non-unique antisense alignments                                                                            
                totalforBlock += mismatchesBlock[x][0][0]+mismatchesBlock[x][0][1]+mismatchesBlock[x][1][0]+mismatchesBlock[x][1][1];
                totalSense += mismatchesBlock[x][0][0] +mismatchesBlock[x][1][0];
                totalAnti += mismatchesBlock[x][0][1] +mismatchesBlock[x][1][1];
                totalUnique += mismatchesBlock[x][0][0]+mismatchesBlock[x][0][1];
                totalUnSense += mismatchesBlock[x][0][0];
                totalUnAnti += mismatchesBlock[x][0][1];
                totalNonUnq += mismatchesBlock[x][1][0]+mismatchesBlock[x][1][1];
                totalNonUnqSense += mismatchesBlock[x][1][0];
                totalNonUnqAnti += mismatchesBlock[x][1][1];
            #summary totals
            totalAlign+=totalforBlock;
            unique+=totalUnique;
           #output totals for all columns                                                      
            output.write(("Totals:").ljust(15)+str(totalforBlock).ljust(20)+str(totalSense).ljust(10)+str(totalAnti).ljust(10));
            output.write(str(totalUnique).ljust(15) +str(totalUnSense).ljust(10) +str(totalUnAnti).ljust(10));
            output.write(str(totalNonUnq).ljust(15) +str(totalNonUnqSense).ljust(10) +str(totalNonUnqAnti).ljust(10)+"\n");
 #check if there is an annotation for this block                                    
            p = cBlockInd * binSize; #get position that this block represents                   
            key = "chr"+chromosomeStr+"."+str(p);
            if annotationsDict.has_key(key):
               #print key, annotationsDict[key];                                               
                output.write(annotationsDict[key]+", " +str(totalforBlock)+ " alignments");
            else:
                output.write("No_Annotation_Found");
            output.write("\n\n\n");
    output.write("TOTAL NUMBER OF ALIGNMENTS: "+str(totalAlign)+"\n");
    output.write("TOTAL UNIQUE ALIGNMENTS: "+str(unique));


#reads all the lines from the different annotation files and returns a list with all of the lines
def readAnnotationFiles(files):
    list =[];
    for x in range(0, len(files), 1):
        print files[x];
        fopen = open(files[x],"r");
        list.extend(fopen.readlines());
        fopen.close();
    return list;

#returns a list with chromosome first, folowed by the starting base pair and ending base pair of the feature
def getChrAndRange(annotation):
    result = [];
    annotationTokens = annotation.split();
    rangeTokens = re.split(r'[=,:,-]', annotationTokens[1]); #splits the string beginning with "range" into tokens
    result.append(rangeTokens[1][3:]);
    result.append(int(rangeTokens[2]));
    result.append(int(rangeTokens[3]));
    return result;

# returns "gene_name$feature$strand"
def getGeneNameAndFeature(annotation, namesDict):
    annotationTokens = annotation.split();
    nameTokens = re.split(r'[_]', annotationTokens[0]);
    geneName = "hello";
    try:
        geneName = namesDict.get(nameTokens[2]+"_"+nameTokens[3]); #finds the well-known name for the gene            
    except:
        print annotation, nameTokens;
    if geneName == None:
        geneName = nameTokens[2];
    concat = geneName+"$"+annotationTokens[5]+"$"+annotationTokens[4];
    return concat;

#returns the bin index of the chromosome if it exists
def getChromNum(chromosome, extraChromList, bins):
    chromNum = -1;
    if chromosome == "X": #to account for X and Y chromosomes                               
        chromNum = 22;
    elif chromosome == "Y":
        chromNum = 23;
    elif chromosome == "M" or len(chromosome) > 2:
        for x in range(0,len(extraChromList),1):
            if extraChromList[x] == chromosome:
                return 23+x+1
    else:
        chromNum = int(chromosome)-1;
    return chromNum;

#iterates through annotations list and creates dictionary of all of the bins contained within that annotation
def annotateBins(annotations, chrBins, extraChromList, binSize, geneDict):
    annotationDict = {};
    for annotation in annotations:
        chrAndRange = getChrAndRange(annotation); #list with chromosome and range values
        chrNum = getChromNum(chrAndRange[0], extraChromList, chrBins);
        if chrNum != -1: #indicates that bin exists for this chromosome
            bin = chrBins[chrNum]; #holds the proper chromosome bin
            start = chrAndRange[1]; #gives the starting bp for this feature 
            numBins = (chrAndRange[2]/binSize)-(chrAndRange[1]/binSize) +1; #endBin's index - startBin's index +1
            for x in range(0, numBins, 1): #need to add extra 1 so loop wont stop prematurely
                p = binSize*(start/binSize)+x*binSize;
                dictID = "chr"+chrAndRange[0]+"."+str(p);
                toAdd =getGeneNameAndFeature(annotation,geneDict);
                if annotationDict.has_key(dictID):
                    existing = annotationDict[dictID];
                    existingTokens = re.split(r'[,]', existing);
                    if toAdd in existingTokens: #weeds out duplicate entries for the same feature
                        continue;
                    else:
                        annotationDict[dictID] =  existing+"," + toAdd;
                else:
                    annotationDict[dictID] = toAdd;
    return annotationDict;


#creates gene names dictionary                                                                    
def createGeneDict(geneNamesFile):
    geneNames = open(geneNamesFile, "r");
    namesList = geneNames.readlines();
    geneNames.close();
    namesDict ={};
    for x in range (1, len(namesList), 1): #start from 1 because first line is header  
        tokens = namesList[x].split();
        namesDict[tokens[0]] = tokens[1];
    return namesDict;
        
#main method
def main(argv):
    inputfile = '';
    outputfile= '';
    binSize = 1000;
    files = [];
    geneNamesFile = '';
    #get the arguments that are passed in
    try:
        opts, args = getopt.getopt(argv,"hi:b:o:f:p:u:t:c:n:g:d:", ["ifile=","collapse"]);
    except getopt.GetoptError:
        print 'Usage: binning_offset_v8.py -i <inputfile> -o <outputfile> -b <binSize> -p <Promoter file> -u <5\'UTR file> -c <CDS file> -t <3\'UTR file> -n <Introns file> -g <gene names file>';
        sys.exit(2);
    for opt, value in opts:
        if opt == '-h':
            print 'Usage: binning_offset_v8.py -i <inputfile> -o <outputfile> -b <binSize> -p <Promoter file> -u <5\'UTR file> -c <CDS file> -t <3\'UTR file> -n <Introns file> -g <gene namesfile>';
            sys.exit();
        #read in all the files
        elif opt in ("-i", "--ifile"):
            inputfile = value;
        elif opt in ("-b"):
            binSize = int(value);
        elif opt in ("-o"):
            outputfile = value;
        elif opt in ("-u"):
            files.append(value);
        elif opt in ("-t"):
            files.append(value);
        elif opt in ("-p"):
            files.append(value);
        elif opt in ("-c"):
            files.append(value);
        elif opt in ("-n"):
            files.append(value);
        elif opt in ("-g"):
            geneNamesFile = value;
    bins= [[] for _ in range(24)]; #create the bins
    extraChromList = [];
    #open file for reading
    input = open(inputfile, "r");
    line = input.readline();
    prevReadName = "string";
    positionMap = {};
    #begins binning process
    while line != "":
        tokens = line.split();
        if tokens[2] == "*": #skips current line and moves to the next alignment if no alignment found
            line = input.readline();
            continue;
        #get chromosome number for alignment
        chromNum = getChromBinNum((tokens[2])[3:], extraChromList, bins);
        startPos = int(tokens[3]);
        
        #continue binning if collapse is false or alignment is new
        mismatches = int((tokens[13])[5:]);
        #get orientation of alignment
        sense = True;
        if tokens[1] == "16":
            sense = False;
        #determine if this alignment comes from a sequencing read that aligns to multiple genomic locations
        multiAlign = False;
        nextLine = input.readline();
        if tokens[0] == prevReadName: #check if read aligned multiple times
            multiAlign = True;
        elif nextLine != "":
            nextLineTokens = nextLine.split();
            if tokens[0] == nextLineTokens[0]:
                multiAlign = True;
        #determine if the genomic-position for the current alignment has been counted already
        pos = str(chromNum)+"."+str(startPos);
        hitNum = positionMap.get(pos);
        #hitNum = 1 means the genomic-position only has had a unique alignment map to it
        #hitNum = 2 means the genomic-position only has had a non-unique alignment map to it
        #hitNum = 3 means the genomic-position has had both a unique and non-unique alignment map to it
        if hitNum != None:
            if multiAlign: #this is a non-unique alignment, so need to check if another non-unique alignment has mapped to it
                if hitNum == 1: #check if position has only seen a unique alignment
                    positionMap[pos] = 3; #update this position since it has now seen both non-unique and unique alignments
                elif hitNum == 2 or hitNum == 3: #check if this position has already see a non-unique alignment; if so, skip to next alignment
                    prevReadName = tokens[0];
                    line = nextLine;
                    continue;
            else: #this is a unique alignment, so need to check if another unique alignment has mapped to it
                if hitNum == 2: #check if position has only seen a non-unique alignment
                    positionMap[pos] = 3; #update this position since it has now seen both unique and non-unique alignments
                elif hitNum == 1 or hitNum == 3: #check if this position has already see a unique alignment; if so, skip to next alignment
                    prevReadName = tokens[0];
                    line = nextLine;
                    continue;
        else: #position has not seen any kind of alignment yet, so set the position code accordingly
            if multiAlign:
                positionMap[pos] = 2;
            else:
                positionMap[pos] = 1;

        #determines which addToBin method to use
        addToBins(chromNum,startPos,bins,binSize, mismatches,multiAlign,sense); #add alignment to proper bin
        prevReadName = tokens[0];
        line = nextLine;

    print "Done binning reads/alignments";
    input.close();
    #get all of the genetic anntotations
    annotations = readAnnotationFiles(files); #creates a huge list of all annotations
    print "Done reading annotation files";
    #create dictionary mapping NCBI IDs for genes to actual gene symbols/names
    namesDict = createGeneDict(geneNamesFile);
    #create annotations dictionary so it is easy to annotate bins
    annotationsDict = annotateBins(annotations, bins, extraChromList, binSize,namesDict);
    print "Done creating annotations dictionary";
    #print the bins with their counts of alignments/insertions and any annotations overlapping with the bins
    output= open(outputfile, "w");
    printBinning(output, bins, binSize, extraChromList,annotationsDict);
    output.close;
    print "Done annotating";

#invoke main method
main(sys.argv[1:]);

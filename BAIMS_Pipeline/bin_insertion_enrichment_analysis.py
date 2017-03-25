"""
Developed by Bhaven Patel, Rajat Rohatgi Lab
This code is free software; you can redistribute it and/or modify it.
Please cite when using or altering this code.
"""

#This script takes the binsAndAnnotations files from multiple datasets and generates determines the anti-sense enrichment bias for bins only annotated as "Intron", insertional enrichment for bins only annotated as "Promoter"; and inactivating insertion enrichment for bins annotated with any genetic feature other than "Promoter". Looks for enrichment only if selected population has at least 1 insertion.

import sys, re, argparse, xlwt
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as compTool


#this class makes it easy to keep track of the information for each bin
class Bin:
    def __init__(self, binTitle,annotation, intronOnly):
        self.binTitle = binTitle; #keeps track of the which chromosome and base pairs this bin is associated with
        self.annotation = annotation; #keeps track of the annotations associated with the bin
        self.intronOnly = intronOnly; #indicator for whether this bin is only annotated with "Intron"
        self.samples = []; #list of the samples that have this bin
        #control stats
        self.controlSense = 0;
        self.controlAnti = 0;
        self.controlUnqSense = 0;
        self.controlUnqAnti = 0;
        self.controlDatasetTotal = 0;
        self.controlDatasetUnqTotal = 0;
        #selected stats
        self.selectedSense = 0;
        self.selectedAnti = 0;
        self.selectedUnqSense = 0;
        self.selectedUnqAnti = 0;
        self.selectedDatasetTotal = 0;
        self.selectedDatasetUnqTotal = 0;
        #p-values
        self.enrichPValue = 0;
        self.adjEnrichPValue = 0;


#this class represents the overall data for an entire sample
class OverallData:
    def __init__(self, sampleName, category, totalSampleInsertions, totalUniqueInsertions):
        self.sampleName = sampleName;
        self.category = category;
        self.totalSampleInsertions = totalSampleInsertions;
        self.totalUniqueInsertions = totalUniqueInsertions;

#this class represents the the specific information for a gene for a sample
class Sample:
    def __init__(self, overallData):
        self.overallData = overallData;
        self.senseInsertions = 0;
        self.antiInsertions = 0;
        self.unqSenInsertions= 0;
        self.unqAntiInsertions = 0;

#This function creates a bin if necessary, populates it with the correct data, and adds it to the specified "binsMap"
def addInsertions(binsMap, binTitle, annotationsLine, intronOnly, overallData, totalsLineTokens, strand):
    bin = None;
    if binsMap.get(binTitle) == None:
        bin = Bin(binTitle, annotationsLine, intronOnly);
        binsMap[binTitle] = bin; #add bin to map of bins
    else:
        bin = binsMap.get(binTitle);
    #create new sample
    sample = Sample(overallData);
    #get sense and anti-sense insertions for this bin in this sample
    if strand: #check which strand the gene is on and assign the insertions accordingly
        sample.senseInsertions = int(totalsLineTokens[2]);
        sample.antiInsertions= int(totalsLineTokens[3]);
        sample.unqSenInsertions = int(totalsLineTokens[5]);
        sample.unqAntiInsertions = int(totalsLineTokens[6]);
    else:
        sample.senseInsertions = int(totalsLineTokens[3]);
        sample.antiInsertions= int(totalsLineTokens[2]);
        sample.unqSenInsertions = int(totalsLineTokens[6]);
        sample.unqAntiInsertions = int(totalsLineTokens[5]);
    bin.samples.append(sample);
    if 'control' in sample.overallData.category: #check if sample is a control dataset
        bin.controlSense += sample.senseInsertions;
        bin.controlAnti += sample.antiInsertions;
        bin.controlUnqSense += sample.unqSenInsertions;
        bin.controlUnqAnti += sample.unqAntiInsertions;
        bin.controlDatasetTotal += sample.overallData.totalSampleInsertions;
        bin.controlDatasetUnqTotal += sample.overallData.totalUniqueInsertions;
    else: #sample is a selected dataset
        bin.selectedSense += sample.senseInsertions;
        bin.selectedAnti += sample.antiInsertions;
        bin.selectedUnqAnti += sample.unqAntiInsertions;
        bin.selectedUnqSense += sample.unqSenInsertions;
        bin.selectedDatasetTotal += sample.overallData.totalSampleInsertions;
        bin.selectedDatasetUnqTotal += sample.overallData.totalUniqueInsertions;


#fills the intronBins and promoterBins maps with the bins and adds the samples
def generateBins(inputLines, intronBins, promoterBins, restBins, overallData):
    for x in range (0,len(inputLines)-3,10):
        binTitle = inputLines[x].strip('\n');
        annotationsLineTokens = inputLines[x+7].split(); #split so we only get annotations and not the alignments
        featureTokens = re.split(r'[$,]', annotationsLineTokens[0]);
        totalsLineTokens = inputLines[x+6].split(); #get line with the totals for the different types of insertions
        #check if this bin covers an intron
        if "Promoter" in featureTokens:
            #skip bin if annotated with any other feature
            if "5\'UTR" not in featureTokens and "CDS" not in featureTokens and "Intron" not in featureTokens:
                #determine the strand the gene is on
                sense = True;
                index = featureTokens.index("Promoter");
                if featureTokens[(index+1)] == "strand=-":
                    sense = False;
                #add bin to "promoterBins"
                intronOnly = False;
                addInsertions(promoterBins, binTitle, annotationsLineTokens[0], intronOnly, overallData, totalsLineTokens, sense);
        #check if this bin covers an intron
        if "Intron" in featureTokens:
            #determine the strand the gene is on
            sense = True;
            index = featureTokens.index("Intron");
            if featureTokens[(index+1)] == "strand=-":
                sense = False;
            intronOnly = False; #initialize this to false
            #add bin to "intronBins" if it isn't annotated with another feature
            if "Promoter" not in featureTokens and "5\'UTR" not in featureTokens and "CDS" not in featureTokens and "3\'UTR" not in featureTokens:
                intronOnly = True;
                addInsertions(intronBins, binTitle, annotationsLineTokens[0], intronOnly, overallData, totalsLineTokens, sense); #add bin to "intronBins"
            #add bin to "restBins"
            addInsertions(restBins, binTitle, annotationsLineTokens[0], intronOnly, overallData, totalsLineTokens, sense);
        elif "5\'UTR" in featureTokens or "CDS" in featureTokens or "3\'UTR" in featureTokens: #check if bin has one of these features if it doesn't have an intron
            #determine the strand the gene is on
            sense = True;
            if "strand=-" in featureTokens:
                sense = False;
            #add bin to "restBins"
            intronOnly = False;
            addInsertions(restBins, binTitle, annotationsLineTokens[0], intronOnly, overallData, totalsLineTokens, sense);



#get the category, whether control or selected, this file represents
def getFileNameAndCategory(fileName):
    fileTokens = re.split(r'[/]', fileName);
    names = [fileTokens[-1]];
    tokens = re.split(r'[_]', fileTokens[-1]);
    if "control" in fileName or "Control" in fileName:
        names.append(tokens[0]+" control");
    else:
        names.append(tokens[0]);
    return names;

#get the total alignments and total unique alignments from binsAndAnnoations file
def getAlignmentTotals(lines):
    alignments = [0]*2;
    tokens = lines[-1].split();
    alignments[1] = int(tokens[3]); #sets unique alignments
    tokens = lines[-2].split();
    alignments[0] = int(tokens[4]); #sets total alignments
    return alignments;


#prints bins annotated with "Intron," which are sorted by enrichment p-value
def printBins2(sheet, intronBinsList):
    #sets up the headers for the sheet
    sheet.write(0,0, "Bin");
    sheet.write(0,1, "Annotation");
    sheet.write(0,2, "p-value");
    sheet.write(0,3, "FDR-corrected p-value");
    numSamples = len(intronBinsList[0].samples);
    col = 4;
    for x in range(0,numSamples):
        sheet.write(0,col, "Sample");
        sheet.write(0,col+1, "Antisense insertions in the bin for the sample");
        sheet.write(0,col+2, "Total insertions mapped in the sample");
        col +=3 ;
    row = 1;
    #print genes
    for bin in intronBinsList:
        if row == 65000: #prevent overstepping excel line limit
            break;
        sheet.write(row, 0, bin.binTitle);
        sheet.write(row, 1, bin.annotation);
        sheet.write(row, 2, bin.enrichPValue);
        sheet.write(row, 3, bin.adjEnrichPValue);
        col = 4;
        samples = bin.samples;
        for sample in samples:
            sheet.write(row, col, sample.overallData.sampleName)
            sheet.write(row, col+1, sample.unqAntiInsertions);
            sheet.write(row, col+2, sample.overallData.totalUniqueInsertions);
            col+=3;
        row += 1;


#prints bins annotated with "Promoter," which are sorted by enrichment p-value
def printBins3(sheet, promoterBinsList):
     #sets up the headers for the sheet
    sheet.write(0,0, "Bin");
    sheet.write(0,1, "Annotation");
    sheet.write(0,2, "p-value");
    sheet.write(0,3, "FDR-corrected p-value");
    numSamples = len(promoterBinsList[0].samples);
    col = 4;
    for x in range(0,numSamples):
        sheet.write(0,col, "Sample");
        sheet.write(0,col+1, "Insertions in the bin for the sample");
        sheet.write(0,col+2, "Total insertions mapped in the sample");
        col +=3 ;
    row = 1;
    #print genes
    for bin in promoterBinsList:
        if row == 65000: #prevent overstepping excel line limit
            break;
        sheet.write(row, 0, bin.binTitle);
        sheet.write(row, 1, bin.annotation);
        sheet.write(row, 2, bin.enrichPValue);
        sheet.write(row, 3, bin.adjEnrichPValue);
        col = 4;
        samples = bin.samples;
        for sample in samples:
            sheet.write(row, col, sample.overallData.sampleName)
            sheet.write(row, col+1, (sample.unqAntiInsertions+sample.unqSenInsertions));
            sheet.write(row, col+2, sample.overallData.totalUniqueInsertions);
            col+=3;
        row += 1;

#prints the rest of the bins with enrichment p-value
def printBins4(sheet, restBinsList):
    #sets up the headers for the sheet
    sheet.write(0,0, "Bin");
    sheet.write(0,1, "Annotation");
    sheet.write(0,2, "p-value");
    sheet.write(0,3, "FDR-corrected p-value");
    numSamples = len(restBinsList[0].samples);
    col = 4;
    for x in range(0, numSamples):
        sheet.write(0,col, "Sample");
        sheet.write(0,col+1, "Sense insertions in the bin for the sample");
        sheet.write(0,col+2, "Antisense insertions in the bin for the sample");
        sheet.write(0,col+3, "Inactivating insertions in the bin for the sample");
        sheet.write(0,col+4, "Total insertions mapped in the sample");
        col +=5 ;
    row = 1;
    #print bins
    for bin in restBinsList:
        if row == 65000: #prevent overstepping excel line limit
            break;
        sheet.write(row, 0, bin.binTitle);
        sheet.write(row, 1, bin.annotation);
        sheet.write(row, 2, bin.enrichPValue);
        sheet.write(row, 3, bin.adjEnrichPValue);
        col = 4;
        samples = bin.samples;
        for sample in samples:
            sheet.write(row, col, sample.overallData.sampleName)
            sheet.write(row, col+1, sample.unqSenInsertions);
            sheet.write(row, col+2, sample.unqAntiInsertions);
            if bin.intronOnly:
            	sheet.write(row, col+3, sample.unqSenInsertions);
            else:
            	sheet.write(row, col+3, (sample.unqSenInsertions+sample.unqAntiInsertions));
            sheet.write(row, col+4, sample.overallData.totalUniqueInsertions);
            col+=5;
        row += 1;

#get the proper bins to analyze based on number of insertions in selected population
def getProperBins(intronBins, promoterBins, restBins):
    #get intron bins that have at least 1 insertion in antisense orientation in selected population
    intronBinsList = []
    for bin in intronBins.values():
        if bin.selectedUnqAnti > 0:
            intronBinsList.append(bin);

    #get promoter bins that have at least 1 insertion in selected population
    promoterBinsList = [];
    for bin in promoterBins.values():
        if (bin.selectedUnqAnti + bin.selectedUnqSense) > 0:
            promoterBinsList.append(bin);

    #get rest bins that have at least 1 insertion in selected population
    restBinsList = [];
    for bin in restBins.values():
    	if bin.intronOnly: #if this is an intron bin, then ensure number of sense insertions is >0
    		if bin.selectedUnqSense > 0:
    			restBinsList.append(bin);
        elif (bin.selectedUnqAnti + bin.selectedUnqSense) > 0: #if bin contains exons, sum insertions and then check if sum is >0
            restBinsList.append(bin);

    return intronBinsList, promoterBinsList, restBinsList;

#calculates the p-value for the bin using fisher test by comparing antisense insertions from control and selected populations for introns or sense+antisense insertions for promoters
def calcPValues(intronBinsList, promoterBinsList, restBinsList):
    #first calculate antisense insertion enrichment p-value for intronic bins
    antisensePVals = []; #list to hold the fisher-test p-values for antisense insertion enrichment in bin for control vs selected datasets
    for bin in intronBinsList:
        #compute p-value for antisense insertion enrichment for bin in control vs selected datasets
        oddsratio, pValue = stats.fisher_exact([[bin.controlUnqAnti, bin.controlDatasetUnqTotal], [bin.selectedUnqAnti, bin.selectedDatasetUnqTotal]], alternative="less");
        bin.enrichPValue = pValue;
        antisensePVals.append(pValue);
    
    adjAntisensePVals = compTool.multipletests(antisensePVals, method='fdr_bh', is_sorted=False)[1];
    for i in range(len(intronBinsList)):
        bin = intronBinsList[i];
        bin.adjEnrichPValue = adjAntisensePVals[i];

    #Calculate p-values for promoter bins
    enrichPVals = []; #list to hold the fisher-test p-values for insertion enrichment in bin for control vs selected datasets
    for bin in promoterBinsList:
        #compute p-value for insertion enrichment for bin in control vs selected datasets
        oddsratio, pValue = stats.fisher_exact([[bin.controlUnqSense + bin.controlUnqAnti, bin.controlDatasetUnqTotal], [bin.selectedUnqSense + bin.selectedUnqAnti, bin.selectedDatasetUnqTotal]], alternative="less");
        bin.enrichPValue = pValue;
        enrichPVals.append(pValue);
    
    adjEnrichPVals = compTool.multipletests(enrichPVals, method='fdr_bh', is_sorted=False)[1];
    for i in range(len(promoterBinsList)):
        bin = promoterBinsList[i];
        bin.adjEnrichPValue = adjEnrichPVals[i];

    #Calculate p-values for rest of the bins
    enrichPVals = []; #list to hold the fisher-test p-values for insertion enrichment in bin in control vs selected datasets
    for bin in restBinsList:
        if bin.intronOnly:
          #compute p-value for sense insertion enrichment for "intron" bin in control vs selected datasets
            oddsratio, pValue = stats.fisher_exact([[bin.controlUnqSense, bin.controlDatasetUnqTotal], [bin.selectedUnqSense, bin.selectedDatasetUnqTotal]], alternative="less");
            bin.enrichPValue = pValue;
            enrichPVals.append(pValue);
        else:
            #compute p-value for total insertion enrichment for bin in control vs selected datasets
            oddsratio, pValue = stats.fisher_exact([[bin.controlUnqSense + bin.controlUnqAnti, bin.controlDatasetUnqTotal], [bin.selectedUnqSense + bin.selectedUnqAnti, bin.selectedDatasetUnqTotal]], alternative="less");
            bin.enrichPValue = pValue;
            enrichPVals.append(pValue);
    adjEnrichPVals = compTool.multipletests(enrichPVals, method='fdr_bh', is_sorted=False)[1];
    for i in range(len(restBinsList)):
        bin = restBinsList[i];
        bin.adjEnrichPValue = adjEnrichPVals[i];



#main parses gene file and does calculations for total alignments per feature
def main(argv):
    #get input  files
    parser = argparse.ArgumentParser();
    parser.add_argument("-i", dest="inputs", help="specifies that the following are input binsAndAnnotations files", nargs="*");
    parser.add_argument("-p", dest="prefix", help="specifies the prefix of the output files");
    args = parser.parse_args();
    inputs = args.inputs;
    prefix = args.prefix;
    intronBins = {}; #map to hold all bins only annotated with Intron
    promoterBins = {}; #map to hold all bins only annotated with Promoter
    restBins = {}; #map to hold all bins overlapping with a genetic feature other than those annotated only as Promoter
    #generate bins objects from bins in each file
    for input in inputs:
        file = open(input , "r");
        inputLines = file.readlines();
        file.close();
        names = getFileNameAndCategory(input); #get category of insertions, whether control or selected
        alignments = getAlignmentTotals(inputLines);
        overallData = OverallData(names[0], names[1], alignments[0], alignments [1]);
        generateBins(inputLines, intronBins, promoterBins, restBins, overallData);
        print "Done generating bins for ", input;
        
    #for each list, retrieve the bins that have at least one insertion in the given category
    intronBinsList, promoterBinsList, restBinsList = getProperBins(intronBins, promoterBins, restBins);
    #calculate p-value for anti-sense insertional enrichment
    calcPValues(intronBinsList, promoterBinsList, restBinsList);
    print "Done calculating p-values for orientation bias and enrichment";
    
    workbook = xlwt.Workbook();
    #sort intron bins by antisense enrichment p-value and print    
    print "There are ", len(intronBinsList), "bins with introns";
    intronBinsList.sort(key = lambda bin: bin.enrichPValue);
    sheet = workbook.add_sheet("Antisense_intronic_bins");
    printBins2(sheet, intronBinsList);

    #sort promoter bins enrichment p-value and print
    print "There are ", len(promoterBinsList), "bins with promoters";
    promoterBinsList.sort(key = lambda bin: bin.enrichPValue);
    sheet = workbook.add_sheet("Promoter_bins");
    printBins3(sheet, promoterBinsList);
    
    #sort rest bins by enrichment p-value
    print "There are ", len(restBinsList), "bins withs with Inactivating insertions";
    restBinsList.sort(key = lambda bin: bin.enrichPValue);
    sheet = workbook.add_sheet("Inactivating_bins");
    printBins4(sheet, restBinsList);
    
    print "Done printing bin insertion enrichment analyses";
    workbook.save(prefix+"_Bin_Analysis.xls");




#main invoked
main(sys.argv[1:]);


"""
Developed by Bhaven Patel, Rajat Rohatgi Lab
This code is free software; you can redistribute it and/or modify it.
Please cite when using or altering this code.
"""

#This script takes the binsAndAnnotations files from a control and selected population and generates a list of genes sorted by significance of insertion enrichment.

import sys, re, argparse, math, xlwt
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as compTool



#keeps track of the statistics relevant to the gene
class Gene:
    def __init__(self, name):
        self.name = name; #keeps track of the which chromosome and base pairs this gene is associated with
        self.samples = []; #list of the samples that have this gene
        self.overlappingGenes = [];
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
        #fisher values
        self.enrichFisherValue = 0;
        self.adjEnrichFisherValue = 0;

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


#fills the genesMap with the genes and adds the samples
def generateGenes(list, genesMap, overallData):
    for x in range (0,len(list)-3,10):
        annotationsLineTokens = list[x+7].split(); #split so we only get annotations and not the alignments
        featureTokens = re.split(r'[$,]', annotationsLineTokens[0]);
        geneAdded = []; #use to keep track of genes that the gene's alignments have been added to
        for y in range(0,len(featureTokens)-1,3):
            #this bin does not overlap no genetic features
            if featureTokens[y] == "No_Annotation_Found":
                continue;
            #case if gene hasnt been recorded yet
            gene = None;
            if genesMap.get(featureTokens[y]) == None:
                gene = Gene(featureTokens[y]);
                genesMap[featureTokens[y]] = gene;
            #gene has been recorded
            else:
                gene = genesMap[featureTokens[y]];
            #add alignments to features data in gene
            totalsLineTokens = list[x+6].split();
            
            #keep track of which genes the insertions in this bin have been added to so we dont add the insertions in this bin to a gene twice
            if featureTokens[y] in geneAdded:
                continue;
            else:
                geneAdded.append(gene.name);

            #determine if the gene has been accounted for in this sample
            samples = gene.samples;
            sample = None;
            for smp in samples:
                if smp.overallData.sampleName == overallData.sampleName:
                    sample = smp;
                    break;
            if sample == None:
                #create new sample
                sample = Sample(overallData);
                gene.samples.append(sample);
            #get sense and anti-sense insertions for this gene in this sample dependng on gene orientation
            if featureTokens[y+2] == "strand=+":
                sample.senseInsertions += int(totalsLineTokens[2]);
                sample.antiInsertions += int(totalsLineTokens[3]);
                sample.unqSenInsertions += int(totalsLineTokens[5]);
                sample.unqAntiInsertions += int(totalsLineTokens[6]);
            else:
                sample.senseInsertions += int(totalsLineTokens[3]);
                sample.antiInsertions += int(totalsLineTokens[2]);
                sample.unqSenInsertions += int(totalsLineTokens[6]);
                sample.unqAntiInsertions += int(totalsLineTokens[5]);
            if 'control' in sample.overallData.category: #add to gene control stats
                if featureTokens[y+2] == "strand=+":
                    gene.controlSense += int(totalsLineTokens[2]);
                    gene.controlAnti += int(totalsLineTokens[3]);
                    gene.controlUnqSense += int(totalsLineTokens[5]);
                    gene.controlUnqAnti += int(totalsLineTokens[6]);
                else:
                    gene.controlSense += int(totalsLineTokens[3]);
                    gene.controlAnti += int(totalsLineTokens[2]);
                    gene.controlUnqSense += int(totalsLineTokens[6]);
                    gene.controlUnqAnti += int(totalsLineTokens[5]);
                if gene.controlDatasetTotal == 0:
                    gene.controlDatasetTotal += sample.overallData.totalSampleInsertions;
                    gene.controlDatasetUnqTotal += sample.overallData.totalUniqueInsertions;
            else: #adds to gene's selected dataset stats
                if featureTokens[y+2] == "strand=+":
                    gene.selectedSense += int(totalsLineTokens[2]);
                    gene.selectedAnti += int(totalsLineTokens[3]);
                    gene.selectedUnqAnti += int(totalsLineTokens[6]);
                    gene.selectedUnqSense += int(totalsLineTokens[5]);
                else:
                    gene.selectedSense += int(totalsLineTokens[3]);
                    gene.selectedAnti += int(totalsLineTokens[2]);
                    gene.selectedUnqAnti += int(totalsLineTokens[5]);
                    gene.selectedUnqSense += int(totalsLineTokens[6]);
                if gene.selectedDatasetTotal == 0:
                    gene.selectedDatasetTotal += sample.overallData.totalSampleInsertions;
                    gene.selectedDatasetUnqTotal += sample.overallData.totalUniqueInsertions;
        #updates the overlapping genes list for each gene
        for geneName in geneAdded:
            gene = genesMap[geneName];
            for name in geneAdded:
                if name in gene.overlappingGenes or name == geneName:
                    continue;
                else:
                    gene.overlappingGenes.append(name);

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

#get the total alignments and total unique alignments from genesAndAnnoations file
def getAlignmentTotals(lines):
    alignments = [0]*2;
    tokens = lines[-1].split();
    alignments[1] = int(tokens[3]); #sets unique alignments
    tokens = lines[-2].split();
    alignments[0] = int(tokens[4]); #sets total alignments
    return alignments;



#prints genes sorted by FDR-corrected p-value for insertion enrichment
def printGenes(sheet, genes):
    #sets up the headers for the sheet
    sheet.write(0,0, "Gene");
    sheet.write(0,1, "p-value");
    sheet.write(0,2, "FDR-corrected p-value");
    numSamples = len(genes[0].samples);
    col = 3;
    for x in range(0,numSamples):
        sheet.write(0,col, "Sample");
        sheet.write(0,col+1, "Insertions in the gene in the sample");
        sheet.write(0,col+2, "Total insertions mapped in the sample");
        col +=3 ;
    row = 1;
    #print genes
    for gene in genes:
        if row == 65000: #prevent overstepping excel line limit
            break;
        sheet.write(row, 0, gene.name);
        sheet.write(row, 1, gene.enrichFisherValue);
        sheet.write(row, 2, gene.adjEnrichFisherValue);
        col = 3;
        samples = gene.samples;
        for sample in samples:
            sheet.write(row, col, sample.overallData.sampleName)
            sheet.write(row, col+1, (sample.unqSenInsertions + sample.unqAntiInsertions));
            sheet.write(row, col+2, sample.overallData.totalUniqueInsertions);
            col+=3;
        row += 1;




#returns a list of Gene objects that have at least one observed insertion in the selected population
def getProperGenes(genesMap):
	genesList = []; #craete list to hold valid Gene objects
	for gene in genesMap.values():
		if (gene.selectedUnqSense + gene.selectedUnqAnti) > 0: #check if total number of insertions in selected population >0
			genesList.append(gene);
	return genesList;



#calculates the p-value for the gene using fisher test by comparing sense/anti-sense values from control and selected populations
def calcPValues(genesList):
    enrichPVals = []; #list to hold the fisher-test p-values for insertional enrichment in gene in selected vs control datasets
    for gene in genesList:
        #calculate enrichment score
        oddsratio, pValue = stats.fisher_exact([[gene.controlUnqSense + gene.controlUnqAnti, gene.controlDatasetUnqTotal], [gene.selectedUnqSense + gene.selectedUnqAnti, gene.selectedDatasetUnqTotal]], alternative="less");
        gene.enrichFisherValue = pValue;
        enrichPVals.append(pValue);
    #get FDR-correct p-values
    adjEnrichPVals = compTool.multipletests(enrichPVals, method='fdr_bh', is_sorted=False)[1];
    #assign each gene its FDR-corrected p-values
    for i in range(len(genesList)):
    	gene = genesList[i];
    	gene.adjEnrichFisherValue = adjEnrichPVals[i];




#main parses gene file and does calculations for total alignments per feature
def main(argv):
    #get input  files
    parser = argparse.ArgumentParser();
    parser.add_argument("-i", dest="inputs", help="specifies that the following are input binsAndAnnotations files", nargs="*");
    parser.add_argument("-p", dest="prefix", help="specifies the prefix of the output files");
    args = parser.parse_args();
    inputs = args.inputs;
    prefix = args.prefix;
    genes = {}; #map to hold all of the genes
    #generate genes from bins in each file
    for input in inputs:
        file = open(input , "r");
        inputLines = file.readlines();
        file.close();
        names = getFileNameAndCategory(input);
        alignments = getAlignmentTotals(inputLines);
        overallData = OverallData(names[0], names[1], alignments[0], alignments [1]);
        generateGenes(inputLines, genes, overallData);
        print "Done generating genes for ", input;

    #get genes with at least one observed insertion in selected population
    genesList = getProperGenes(genes);
    calcPValues(genesList);
    print "Done calculating p-values for gene insertion enrichment";

    workbook = xlwt.Workbook();
    #sort by enrichFisherValue
    genesList.sort(key = lambda gene: gene.adjEnrichFisherValue);
    sheet = workbook.add_sheet("Gene_Insert_Enrich");
    printGenes(sheet, genesList);
    workbook.save(prefix+"_Gene_Comparison.xls");

    print "Done printing gene insertion enrichment p-values";




#main invoked
main(sys.argv[1:]);


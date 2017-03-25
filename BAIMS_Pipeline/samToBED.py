"""
Developed by Bhaven Patel, Rajat Rohatgi Lab
Modified from code written by Aaron Quinlan on 2009-08-27.
This code is free software; you can redistribute it and/or modify it.
Please cite when using or altering this code.
"""

#This script converts the uniquely aligned reads from the .sam file to a .bed file format. Only one alignment per genome position will be written to the .bed file.

import sys
import getopt
import re

#load .sam file and convert uniquely aligned reads to .bed file format; only one alignment per position in the genome will be written to the .bed file
def processSAM(file,outputfile):
    f = open(file,'r');
    output= open(outputfile, 'w');
    prevReadName="";
    positionMap={};
    line = f.readline();
    while line != "":
        tokens = line.split();
        if tokens[2] == "*": #skips current line and moves to the next alignment if no alignment found
        	line = f.readline();
        	continue;
        #get chromosome
        chromNum = tokens[2];
        #get sequence
        sequence = tokens[9];
        #get start position and end position
        startPos = tokens[3];
        endPos = str(int(startPos)+len(sequence)-1);
        if int(tokens[1]) == 16:
            sense = "-";
        else:
            sense ="+";
        #check if read is unique
        multiAlign = False;
        nextLine = f.readline();
        if tokens[0] == prevReadName: #check if alignment comes from the same sequencing read as the previous alignment
            multiAlign = True;
        elif nextLine != "":
            nextLineTokens = nextLine.split();
            if tokens[0] == nextLineTokens[0]: #check if alignment comes from the same sequencing read as the next alignment
                multiAlign = True;
        #found unique read
        if not multiAlign:
            #check if position has been hit
            pos = chromNum+"."+startPos;
            hitNum = positionMap.get(pos);
            if hitNum == None:
                #put position in positionMap
                positionMap[pos]=1;
                #write position in output file
                entry = chromNum + "\t" + startPos + "\t" + endPos + "\t"+sequence+"\t"+"0"+"\t"+sense+"\n";
                output.write(entry);
        #update prevReadName
        prevReadName=tokens[0];
        line = nextLine;

    f.close();
    output.close();


#main method
def main(argv=None):
    if argv is None:
        argv = sys.argv;
    opts, args = getopt.getopt(argv[1:], "hi:o:", ["help"]);

    # option processing
    samFile = "";
    outputfile = "";
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i <input a .SAM file> -o <any output file name>";
        if option in ("-i"):
            samFile = value;
        if option in ("-o"):
            outputfile = value;
                    
    # make a BED file of the SAM file.
    print "making BED file";
    processSAM(samFile, outputfile);
    print "done making BED file";

if __name__ == "__main__":
    sys.exit(main());
#!/bin/sh

#Bhaven Patel
#Wrapper script for running haploid screen pipeline. This pipeline uses BAIMS to analyze a control and selected population for both genes that are enriched for GT insertions and bins that are enriched for GT insertions hypothesized not to affect gene function.

#OPTIONS
# -c <file containing reads for the control/unselected population>
# -s <file containing reads for the sorted/selected population>
# --binSize <followed by the size of the bins to hold the alignments> This option only works with a fastq file when the bins are created.
# -r <specifies that binning has already been run on the control and selected files>
# -n <name for output files>

#NEED TO COMPLETE THIS WITH THE PATH TO THE FOLDER THAT HOLDS THE BOWTIE-1.0.1 EXECUTABLE; DO NOT INCLUDE "/" FOLLOWING THE NAME OF THE FOLDER
bowtie_Path="path_to_Bowtie_folder";
echo "Path to the Bowtie folder is ${bowtie_Path}";
#NEED TO COMPLETE THIS WITH THE PATH TO THE BOWTIE INDEX TO USE FOR ALIGNMENT, AS SPECIFIED IN THE BOWTIE MANUAL
bowtie_index_Path="path_to_Bowtie_index_for_human_genome_GRCh38";
echo "Path to the Bowtie index needed for alignment is ${bowtie_index_Path}";

#NEED TO COMPLETE THIS WITH THE PATH TO THE "BAIMS_Pipeline" FOLDER; DO NOT INCLUDE "/" FOLLOWING THE NAME OF THE FOLDER
pipeline_Path="path_to_BAIMS_Pipeline_folder";
echo "Path to the BAIMS_pipeline folder is ${pipeline_Path}";

control_file=;
selected_file=;
file_prefix="output";
binSize=1000;
alreadyRun=false;

#get all options                                                                       
while [ $# -gt 0 ]
do
    case "$1" in
        -c) control_file="$2"; shift;;
        -s) selected_file="$2"; shift;;
        -r) alreadyRun=true;;
        --binSize) binSize="$2"; shift;;
		-n) file_prefix="$2"; shift;;
        *) echo "$1 is not a valid option.";
                exit 1;;
    esac
    shift
done

#check that bin size is an integer
re='^[0-9]+$';
if ! [[ $binSize =~ $re ]]
then
	echo "Need to provide an integer for the bin size";
	exit 1;
fi


#check that control and selected fastq files have been given
if [[ -z "$control_file" || -z "$selected_file" ]]
then
	echo "Need to give the path to the fastq files for both the control and selected populations";
	exit 1;
fi

#loop through and run bowtie and binning on control and selected files if they have not been processed before
if [ $alreadyRun == false ]
then
	i=0;
	until [ $i -gt 1 ]
	do
		#do Bowtie Alignment
		fastq="";
		prefix="";
		if [ $i -eq 0 ]
		then
			fastq="$control_file";
			prefix="${file_prefix}_control";
		else
			fastq="$selected_file";
			prefix="${file_prefix}_selected";
		fi

		outputfile="$prefix.sam";
		echo "Beginning bowtie alignment for $fastq";
		${bowtie_Path}/bowtie -S --sam-nohead -p 4 -v 3 -k 5 --best ${bowtie_index_Path} $fastq $outputfile

		#do initial analysis of the number of reads aligned, the number of mismatches for aligned reads, and to which chromosomes reads align
		inputfile=$outputfile;
		outputfile="${prefix}_alignment_stats";
		python ${pipeline_Path}/initial_sam_analysis.py -i $inputfile -o $outputfile
		alignmentsFile=$inputfile;

		#create BED file
		bedOutput="$prefix.bed";
    	python ${pipeline_Path}/samToBED.py -i $alignmentsFile -o $bedOutput;

    	#create bins file
   		binFile="${prefix}_binsAndAnnotations";
		python ${pipeline_Path}/binning.py -i $alignmentsFile -o $binFile -b $binSize -p ${pipeline_Path}/GRCh38_annotations/GRCh38_Promoters -u ${pipeline_Path}/GRCh38_annotations/GRCh38_5\'UTR_Exons -t ${pipeline_Path}/GRCh38_annotations/GRCh38_3\'UTR_Exons -c ${pipeline_Path}/GRCh38_annotations/GRCh38_CDS_Exons -n ${pipeline_Path}/GRCh38_annotations/GRCh38_Introns -g ${pipeline_Path}/GRCh38_annotations/GRCh38_geneNames
		
		#update index
		i=`expr $i + 1`;
	done
fi

#create input line for bin and gene analysis
inputs="${file_prefix}_control_binsAndAnnotations ${file_prefix}_selected_binsAndAnnotations";
echo $inputs;

echo "WORKING ON BIN ANALYSES";
python ${pipeline_Path}/bin_insertion_enrichment_analysis.py -i $inputs -p "$file_prefix";
echo "FINISHED BIN ANALYSES";


echo "WORKING ON GENE ANALYSIS";
python ${pipeline_Path}/gene_insertion_enrichment_analysis.py -i $inputs -p "$file_prefix" ;
echo "FINISHED GENE ANALYSIS";

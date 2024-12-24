#!/bin/bash

###################################################################################################################################
# This script has been used to search for ARM, VF and MGEs on the metagenomics reads with the program abricate more info:
<https://github.com/tseemann/abricate>

# It is important to note that it has been designed for a specific working directory. Therefore, the reproduction of the results
# will require small modifications of the script or the adaptation of your working directory.

# Created on Thursday September 1st 2024.

# @author: PhD C. Julio Cesar Ortega Camabara - School of Biomedical Sciences, UWL.

# Version: 1
####################################################################################################################################

# ENV VARIABLES
# source /home/julio92-c/anaconda3/etc/profile.d/conda.sh

# Create the following list:
# -DIR contains the dataset information. The structure of the information contained in DIR is: (Name). For example: (Escherichia_coli).
DIR=(A60B A61)

DBs=(ncbi card plasmidfinder vfdb)

echo -e "\n\n"---------------------------STARTING----------------------------

for ((i=0; i<=${#DIR[@]}-1; i++)); do


        # Check the existence of directory.
        FOLDER=${DIR[$i]}
        if [ -d $FOLDER ]; then
                echo -e "\n\n"DIR: $FOLDER
                echo
                cd $FOLDER
        else
                echo $FOLDER "doesn't exist"
        echo -e "\n\n"Creating dir $FOLDER ...
        echo $FOLDER was succefully created
        mkdir $FOLDER
                echo -------------------------------------------------------------------------------------------------------
                continue
        fi

	# Check the existence of dataset file.
         # The dataset file name has always the same structure: dataset name like (dir_name).fasta
        #FILE=${DIR[$i]}_genomic_reference.gb
        FILE1=${DIR[$i]}_se_db1.kraken2.classified.fastq.gz
	FILE2=${DIR[$i]}_kraken2ID_PFP.fastq.gz
	FILE3=${DIR[$i]}_kraken2ID_PFP.fasta

	echo "Merging Kraken tax ID with read ID..."
        gunzip -c $FILE1 | sed '/^@/s/ /_/g' | gzip > $FILE2

	source activate abricate

	echo "Converting fastq.gz to fasta format..."
	any2fasta $FILE2 > $FILE3

	echo "Removing temporary files..."
	rm $FILE2


	# AMR profiling
	for ((j=0; j<=${#DBs[@]}-1; j++)); do

	db=${DBs[$j]}


        if [ -s $FILE3 ]; then

                # Run stats on Genomic Dataset
		echo -e "\n\n"DB: $db
                echo -e "\n\n"FILE: $FILE3
		
		echo "Starting ARM annotation..."
		abricate --version
		abricate $FILE3 --db $db --minid 70 --mincov 70 --csv > ${DBs[$j]}.csv

        else

		 echo $FILE2 "doesn't exist"

         fi
	done
	
	echo "Closing the env..."
	conda deactivate
	echo "Generating the summary report..."
	cat *.csv > $FOLDER_summary_report.csv
        echo -----------------------------------------------------------------
        cd ../
done

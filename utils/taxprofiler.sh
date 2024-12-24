#!/bin/bash -l

# Batch script to run an OpenMP threaded job under SGE.

#$ -m be
#$ -M 32104617@student.uwl.ac.uk

# Request five hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=8:30:0

# Request 2 gigabyte of RAM for each core/thread
# (must be an integer followed by M, G, or T)
#$ -l mem=2G

# Request 100 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=100G

# Set the name of the job.
#$ -N taxprofiler

# Request 16 cores.
#$ -pe smp 16

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID
#$ -wd /home/rekgort/Scratch/taxprofiler-tutorial

# 8. Run the application.
nextflow run nf-core/taxprofiler -r 1.1.8 -profile singularity \
--input ./samplesheet.csv \
--databases ./database.csv \
--outdir ./oscar_results \
--perform_longread_qc \
--perform_longread_hostremoval --hostremoval_reference ./GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
--perform_runmerging --save_runmerged_reads \
--run_kraken2 \
--kraken2_save_readclassifications --kraken2_save_reads \
--run_profile_standardisation \
--run_krona \
-resume
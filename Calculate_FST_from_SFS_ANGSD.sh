#!/bin/bash
#SBATCH --clusters=all
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --array=1-26:1
#SBATCH --time=7-00:00:00 
#SBATCH --job-name=CalcFST
#SBATCH --partition=long
#SBATCH --output=CalcFST_%A_%a.log
#SBATCH --error=CalcFST_%A_%a.error


###############################################################################################################
# CALCULATE FST FROM SFS USING ANGSD AND REALSFS
# A.Sendell-Price, Nov 2021
###############################################################################################################

#Notes: 
#The following script can be used to calculate genome-wide and window based estimates of FST between
#a pair of populations based on either folded or un-folded Site Frequency Spectra SFS using the programmes
#angsd (version 0.935) and realSFS.

#Prior to submitting the script you will need to update the number of jobs in the sbatch array command at the
#top of the script matches the number of chromosomes for your species. For example, for 26 chromosomes this
#line needs to read as follows:
#SBATCH --array=1-26:1
#whereas for 33 chromosomes the line will need to read:
#SBATCH --array=1-33:1

#The input files in the following code chunk will need to be updated, providing the following information:
#1) required window size and step size in bp (currently set to 100000)
#2) the absolute path to a chromosome list, specifying the chromosome names (1 per line) exactly as they appear
#in the reference assembly
#3) the absolute path to a list of bam files for each of the two populations you want to compare
#4) the absolute path to the reference assembly for your species
#5) either the absolute path to the a consensus sequence for the outgroup (if you want to use the un-folded SFS)
#or set ANC to REF (if you want to use the folded SFS)
#6) A name for your output directory e.g. popX_vs_popY


#############################################################################################################
# DEFINE SETTING FOR ANALYSIS (CHANGE TO MATCH YOUR FILES)
#############################################################################################################

#Define window and step size to use
WINDOW_SIZE=100000
STEP_SIZE=100000

#Set path to list of chromosome names (needs to match chrom names exactly as they appear in reference assembly)
CHROM_LIST=/path/to/chrom/list

#Set path to list of sample bams for each population, this will need to be the absolute path
POP1_BAM_LIST=/path/to/pop1/bam/list
POP2_BAM_LIST=/path/to/pop2/bam/list

#Set path to the Reference genome (ends *fasta.gz)
REF=/path/to/reference/assembly

#Either set path to ancestral consensus sequence (to calcualte derived SFS aka unfolded SFS) or if you dont have this set ancestral to referene (the folded SFS will be calculated instead)
#ANC=/path/to/ancestral/consensus/sequence #<- this version for unfolded SFS
ANC=$REF #<- this version for folded SFS

#Set output directory name (something sensible like PopX_vs_PopY)
DIRECTORY=name_of_directory


#############################################################################################################
# LEAVE FOLLOWING COMMANDS AS THEY ARE
#############################################################################################################

#Load required module
ml angsd/0.935-GCC-10.2.0

#Get chromosome name using slurm job array ID
CHROM=$(cat $CHROM_LIST | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#Create output directory
mkdir $DIRECTORY
cd $DIRECTORY

#Calculate FST based on folded SFS if REFERENCE and ANCESTRAL match . . . 
if [ $REF = $ANC ]
then
    #Make a directory for the chromosome and move into it
    mkdir $CHROM
    cd $CHROM

    #Calculate the allele frequency likelihoods for pop1 and pop2 using angsd
    angsd -bam $POP1_BAM_LIST -ref $REF -anc $ANC -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 20 -minQ 20 -gl 1 -doSaf 1 -r $CHROM -out pop1.${CHROM}
    angsd -bam $POP2_BAM_LIST -ref $REF -anc $ANC -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 20 -minQ 20 -gl 1 -doSaf 1 -r $CHROM -out pop2.${CHROM}

    #Calculate the 2D site frequency spectrum
    #-fold 1: 			tells realSFS to output the "folded" site frequency spectrum
    #-m 0: tells		tells realSFS to use the standard EM algorithm instead of the accelerated EM algorithm
    #-maxIter 50000: 	sets the maximum number of EM iterations to 50,000 instead of the default 100 (to provide sufficient iterations for covergence)
    #-tole 1e-6: 		when the difference in successive likelihood values in the EM algorithm gets below this value the optimisation will stop.
    realSFS pop1.${CHROM}.saf.idx pop2.${CHROM}.saf.idx -fold 1 -m 0 -maxIter 50000 -tole 1e-4 > pop1.pop2.${CHROM}.folded.sfs

    #Calculate FST on a per site basis
    realSFS fst index pop1.${CHROM}.saf.idx pop2.${CHROM}.saf.idx -sfs pop1.pop2.${CHROM}.folded.sfs -fstout pop1.pop2.${CHROM}

    #Get the global FST estimate for chromosome and save output to file
    realSFS fst stats pop1.pop2.${CHROM}.fst.idx > ${CHROM}.global.fst

    #Estimate FST in a sliding window
    #(-type 2 sets the left most position of the first window to 1 i.e. the first bp of the chromosome)
    realSFS fst stats2 pop1.pop2.${CHROM}.fst.idx -win $WINDOW_SIZE -step $STEP_SIZE -type 2 > pop1.pop2.${CHROM}.fst.size${WINDOW_SIZE}_step${STEP_SIZE}
    cd ../
fi


#Calculate FST based on unfolded SFS if REFERENCE and ANCESTRAL don't match . . . 
if [ $REF != $ANC ]
then
    #Make a directory for the chromosome and move into it
    mkdir $CHROM
    cd $CHROM

    #Calculate the allele frequency likelihoods for pop1 and pop2 using angsd
    angsd -bam $POP1_BAM_LIST -ref $REF -anc $ANC -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 20 -minQ 20 -gl 1 -doSaf 1 -r $CHROM -out pop1.${CHROM}
    angsd -bam $POP2_BAM_LIST -ref $REF -anc $ANC -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -minMapQ 20 -minQ 20 -gl 1 -doSaf 1 -r $CHROM -out pop2.${CHROM}

    #Calculate the 2D site frequency spectrum
    #-fold 0: 			tells realSFS to output the "unfolded" site frequency spectrum
    #-m 0: tells		tells realSFS to use the standard EM algorithm instead of the accelerated EM algorithm
    #-maxIter 50000: 	sets the maximum number of EM iterations to 50,000 instead of the default 100 (to provide sufficient iterations for covergence)
    #-tole 1e-6: 		when the difference in successive likelihood values in the EM algorithm gets below this value the optimisation will stop.
    realSFS pop1.${CHROM}.saf.idx pop2.${CHROM}.saf.idx -fold 0 -m 0 -maxIter 50000 -tole 1e-4 > pop1.pop2.${CHROM}.unfolded.sfs

    #Calculate FST on a per site basis
    realSFS fst index pop1.${CHROM}.saf.idx pop2.${CHROM}.saf.idx -sfs pop1.pop2.${CHROM}.unfolded.sfs -fstout pop1.pop2.${CHROM}

    #Get the global FST estimate for chromosome and save output to file
    realSFS fst stats pop1.pop2.${CHROM}.fst.idx > ${CHROM}.global.fst

    #Estimate FST in a sliding window
    #(-type 2 sets the left most position of the first window to 1 i.e. the first bp of the chromosome)
    #-type 2: 
    realSFS fst stats2 pop1.pop2.${CHROM}.fst.idx -win $WINDOW_SIZE -step $STEP_SIZE -type 2 > pop1.pop2.${CHROM}.fst.size${WINDOW_SIZE}_step${STEP_SIZE}
    cd ../
fi

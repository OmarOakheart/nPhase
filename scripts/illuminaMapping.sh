#!/bin/bash

#Initializing
strainName=$1 #PDB1_S1
R1=$2 #${rawIlluminaReadFolder}/PDB1_S1_R1_001.fastq.gz
R2=$3 #${rawIlluminaReadFolder}/PDB1_S1_R2_001.fastq.gz
reference=$4 #$reference
outputFolder=$5 #$mappedReadFolder

#Mapping
samFile=${outputFolder}/${strainName}.sam
RGline="@RG\tID:ID_${strainName}\tLB:LB_${strainName}\tPL:ILLUMINA\tPU:PU_${strainName}\tSM:SM_${strainName}"
bwa mem -M -R $RGline $reference $R1 $R2 > $samFile

#Organizing mapped reads
sortedBamFile=${outputFolder}/${strainName}.sorted.bam

echo "Removing quality 0 (multimapped) reads, turning to bam and sorting it"
samtools view -bT $reference -q 1 $samFile | samtools sort - > $sortedBamFile
samtools index $sortedBamFile

#GATK cleaning
MDsortedBamFile=${outputFolder}/${strainName}.sorted.MD.bam
gatk MarkDuplicates --REMOVE_DUPLICATES true -I $sortedBamFile -O $MDsortedBamFile -M /dev/null

#Finalizing
finalBamFile=${outputFolder}/${strainName}.final.bam
samtools sort $MDsortedBamFile -o $finalBamFile
samtools index $finalBamFile

#Cleaning up
rm $samFile $sortedBamFile $MDsortedBamFile


#!/bin/bash

#Initializing
strainName=$1
choppedReads=$2
reference=$3
outputFolder=$4
flag=$5
threads=$6

#Map with ngmlr
ngmlr -t $threads -x ont -r $reference -q $choppedReads -o ${outputFolder}/${strainName}.sam

#Convert sam to bam
samtools view -bT $reference -F $flag ${outputFolder}/${strainName}.sam | samtools sort - > ${outputFolder}/${strainName}.sorted.bam
samtools index ${outputFolder}/${strainName}.sorted.bam
samtools view ${outputFolder}/${strainName}.sorted.bam > ${outputFolder}/${strainName}.sorted.sam

#Cleaning up
rm ${outputFolder}/${strainName}.sam


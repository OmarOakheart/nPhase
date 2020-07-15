#!/bin/bash

#Required user-defined parameters
strainName=$1
illuminaFile_R1=$2
illuminaFile_R2=$3
minIONFile=$4
reference=$5 #/home/oabousaada/ReferenceSequence/Sace/R64_nucl.fasta
outputFolder=$6 #/home/oabousaada/HPhasing/phenovarPolyploids
threads=$7

#These are default parameters, they should only be modified if overwritten
maxID=0.05
minOvl=0.1
minLen=0
minSim=0.01

#Need a way to make sure nPhase can always access its scripts and the correct python version
scriptFolder=/home/oabousaada/HPhasing/Scripts
python=/home/oabousaada/Soft/miniconda3/envs/haplotyping/bin/python3

#Create folder structure
basePath=${outputFolder}/${strainName}
mkdir -p ${basePath}/Mapped/Illumina
mkdir -p ${basePath}/VariantCalls/Illumina
mkdir -p ${basePath}/Mapped/MinION
mkdir -p ${basePath}/VariantCalls/MinION
mkdir -p ${basePath}/Overlaps
mkdir -p ${basePath}/Phased
mkdir -p ${basePath}/Logs

#Links to scripts
minIONMappingScript=${scriptFolder}/minIONMapping.sh
illuminaMappingScript=${scriptFolder}/illuminaMapping.sh
getIlluminaPositionScript=${scriptFolder}/getIlluminaPositions.py
SNPminIONAssignmentScript=${scriptFolder}/assignMinIONToSNPs4_2.py
minIONValidationScript=${scriptFolder}/MinIONValidation.py
haplotypeResolverScript=/home/oabousaada/HPhasing/Scripts/nPhase.py

##########################
#Pre-process minION reads#
##########################

#Remove adapters with porechop
choppedMinION=${basePath}/rawReads/MinION/porechop_${minIONFile}
porechop -i $minIONFile -o $choppedMinION --threads $threads

#Map porechopped reads to reference
mappedMinIONFolder=${basePath}/Mapped/MinION
minIONFlag=260 #This flag allows split reads
sh $minIONMappingScript $strainName $choppedMinION $reference $mappedMinIONFolder $minIONFlag $threads

############################
#Pre-process Illumina reads#
############################

#Map Illumina reads to reference
mappedIlluminaFolder=${basePath}/Mapped/Illumina
sh $illuminaMappingScript $strainName $illuminaFile_R1 $illuminaFile_R2 $reference $mappedIlluminaFolder

#Variant call Illumina reads and select SNPs only
calledIlluminaFolder=${basePath}/VariantCalls/Illumina
estimatedPloidy=2 #This argument is required for GATK's variant calling and we can expect most SNPs to only have two alleles anyway
gatk HaplotypeCaller -R $reference -ploidy $estimatedPloidy -I ${mappedIlluminaFolder}/${strainName}.final.bam -O ${calledIlluminaFolder}/${strainName}.vcf
gatk SelectVariants -R $reference --variant ${calledIlluminaFolder}/${strainName}.vcf -O ${calledIlluminaFolder}/${strainName}.SNPs.vcf --select-type-to-include SNP

#Extract heterozygous positions from VCF file
cat ${calledIlluminaFolder}/${strainName}.SNPs.vcf | grep "#" >${calledIlluminaFolder}/${strainName}.hetSNPs.vcf
cat ${calledIlluminaFolder}/${strainName}.SNPs.vcf | grep -v "#" | sed 's/;/\t/g' | grep -v "AF=1.00" | awk -vOFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8}' >>${calledIlluminaFolder}/${strainName}.hetSNPs.vcf
vcfFile=${calledIlluminaFolder}/${strainName}.hetSNPs.vcf
illuminaPositionsOutputFile=${mappedIlluminaFolder}/${strainName}.hetSNPs.positions.tsv
$python $getIlluminaPositionScript $vcfFile $illuminaPositionsOutputFile
illuminaSNPsBed=${mappedIlluminaFolder}/${strainName}.hetSNPs.bed
cat $illuminaPositionsOutputFile | sed 's/:/\t/g' | awk -vOFS="\t" '{print $1,$2-1,$2}' >$illuminaSNPsBed

#Reduce long reads to their heterozygous SNPs
cleanMinIONSamFile=${basePath}/Mapped/MinION/${strainName}.sorted.sam
minQ=0.01   # Currently
minMQ=0     # Not
minAln=0.5  # Used
minIONPositionNTFile=${basePath}/VariantCalls/MinION/${strainName}.hetPositions.SNPxMinION.tsv
$python $SNPminIONAssignmentScript $cleanMinIONSamFile $illuminaSNPsBed $reference $minQ $minMQ $minAln $minIONPositionNTFile

#Only keep the longest reads and get rid of duplicate heterozygous SNP profiles, compensate by keeping context coverage information.
minCov=0	 # Currently
minRatio=0	 # Not
minTrioCov=0 # Used
validatedSNPAssignmentsFile=${basePath}/VariantCalls/MinION/${strainName}.hetPositions.SNPxMinION.validated.tsv
contextDepthsFile=${basePath}/Overlaps/${strainName}.contextDepths.tsv
$python $minIONValidationScript $minIONPositionNTFile $minCov $minRatio $minTrioCov $validatedSNPAssignmentsFile $contextDepthsFile

#Run everything through nPhase
phasedFolder=${basePath}/Phased
$python $haplotypeResolverScript $validatedSNPAssignmentsFile $strainName $contextDepthsFile $phasedFolder $reference $minSim $minOvl $minLen $maxID

#Generate fastQ files from haplotig clusters

#Need to write that script, the file you want is 
haplotigReadNameFile=${phasedFolder}/${strainName}_${minOvl}_${minSim}_${maxID}_clusterReadNames.tsv


#Code that tries to generate graphs for datavis?



#Code that echoes the different files and how to find them




import sys

minIONFilePath=sys.argv[1]
minCov=int(sys.argv[2])
minRatio=float(sys.argv[3])
minTrioCov=int(sys.argv[4])
validatedSNPAssignmentsFile=sys.argv[5]
contextDepthsFile=sys.argv[6]

#minIONFilePath="/home/oabousaada/HPhasing/WAxNA/VariantCalls/MinION/WAxNA.hetPositions.MinIONxSNPs.tsv"
#minCov=0
#minRatio=0
#minTrioCov=0
#validatedSNPAssignmentsFile="/home/oabousaada/HPhasing/WAxNA/VariantCalls/MinION/WAxNA.hetPositions.MinIONxSNPs.validated.tsv"
#contextDepthsFile="/home/oabousaada/HPhasing/WAxNA/Overlaps/WAxNA.contextDepths.tsv"

minIONFile=open(minIONFilePath,"r")

minIONReads=[]

for read in minIONFile:
	minIONReads.append(read.strip("\n").split("\t"))

minIONFile.close()

contextDepths={}
distance=2

for line in minIONReads:
	SNPs=line[1:]
	for SNP in SNPs:
		i=SNPs.index(SNP)
		itemList=[]
		for j in range(i-distance,i+distance+1):
			if j>=0 and j<len(SNPs):
				itemList.append(SNPs[j])
		if SNP not in contextDepths:
			contextDepths[SNP]={frozenset(itemList):0}
		else:
			contextDepths[SNP][frozenset(itemList)]=0

for line in minIONReads:
	SNPs=line[1:]
	for SNP in SNPs:
		i=SNPs.index(SNP)
		itemList=[]
		for j in range(i-distance,i+distance+1):
			if j>=0 and j<len(SNPs):
				itemList.append(SNPs[j])
		contextDepths[SNP][frozenset(itemList)]+=1

contextText=""

for SNP, contexts in contextDepths.items():
	contextTextLine=SNP+"\t"
	for contextSNPs, quantity in contexts.items():
		contextTextLine+=" ".join(contextSNPs)+"\t"+str(quantity)+"\t"
	contextTextLine=contextTextLine.strip("\t")
	contextText+=contextTextLine+"\n"

contextFile=open(contextDepthsFile,"w")
contextFile.write(contextText)
contextFile.close()


minIONSets={}
minIONTuples={}


###########################################################################################
#You can't just do this step before you remove reads based on coverage, that's ridiculous.#
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
###########################################################################################

#Only keeping MinION reads that have unique combinations of SNPs 
uniqueMinIONReadSets=set()
uniqueMinIONReadIDs=set()

curLen=len(uniqueMinIONReadSets)

for read in minIONReads:
	readSet=frozenset(read[1:])
	altruisticSet=[]
	for SNP in readSet:
		altruisticSet.append(SNP)
	readSet=frozenset(altruisticSet)
	uniqueMinIONReadSets.add(readSet)
	if len(uniqueMinIONReadSets)>curLen:
		curLen+=1
		minIONSets[read[0]]=readSet
		minIONTuples[read[0]]=tuple(readSet)
		uniqueMinIONReadIDs.add(read[0])

###########################################################################################
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#You can't just do this step before you remove reads based on coverage, that's ridiculous.#
###########################################################################################

SNPAssignments=minIONSets

#Printing it all to a file
SNPAssignmentsText=""

for longRead, assignedSNPs in SNPAssignments.items():
	if assignedSNPs!=set({''}):
		SNPAssignmentsText+=longRead+"\t"+"\t".join(assignedSNPs)+"\n"

SNPAssignmentFile=open(validatedSNPAssignmentsFile,"w")
SNPAssignmentFile.write(SNPAssignmentsText)
SNPAssignmentFile.close()

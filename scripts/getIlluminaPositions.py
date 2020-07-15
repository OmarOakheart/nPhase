import sys

#Retrieving arguments
hetVCFName=sys.argv[1]
illuminaPosFileName=sys.argv[2]

#hetVCFName="/home/oabousaada/HPhasing/WAxNA/VariantCalls/Illumina/WAxNA.hetSNPs.vcf"
#illuminaPosFileName="/home/oabousaada/HPhasing/WAxNA/Mapped/Illumina/WAxNA.hetSNPs.positions.tsv"

hetVCFFile=open(hetVCFName,"r")

hetVCFPositions=set()

for line in hetVCFFile:
	line=line.strip("\n").split("\t")
	if "#" not in line[0]:
		chr=line[0]
		pos=line[1]
		ref=line[3]
		alts=line[4].split(",")
		for alt in alts:
			if len(alt)==1 and len(ref)==1: #SNP or del of 1 base ("*")
				hetVCFPositions.add(chr+":"+pos)
			elif len(alt)>len(ref): #INS
				#hetVCFPositions.add(chr+":"+pos)
				pass
			elif len(alt)<len(ref): #DEL
				#for curPos in range(int(line[1])+1,int(line[1])+len(ref)):
				#	hetVCFPositions.add(chr+":"+str(curPos))
				pass

#Writing the Illumina het positions to a file
illuminaSNPText=""

for position in hetVCFPositions:
	illuminaSNPText+=position+"\n"

illuminaSNPFile=open(illuminaPosFileName,"w")

illuminaSNPFile.write(illuminaSNPText)

illuminaSNPFile.close()

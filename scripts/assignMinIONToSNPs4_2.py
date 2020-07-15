import re
import sys

###########
#Resources#
###########

#Cigar string:
#http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
#http://zenfractal.com/2013/06/19/playing-with-matches/
#https://davetang.org/wiki/tiki-index.php?page=SAM

#QScores
#https://drive5.com/usearch/manual/quality_score.html

#To consider:
#Currently, softclips, insertions and padding worsen the alignment score, you might want to change that

cigarPartRE='(\d+)([DHIMNPSX=])'

samPath=sys.argv[1]
bedPath=sys.argv[2]
referencePath=sys.argv[3]
minQ=float(sys.argv[4])
minMQ=int(sys.argv[5])
minAln=float(sys.argv[6])
outputPath=sys.argv[7]

#samPath="/home/oabousaada/HPhasing/WAxNA/Mapped/MinION/WAxNA.sorted.sam"
#bedPath="/home/oabousaada/HPhasing/WAxNA/Mapped/Illumina/WAxNA.hetSNPs.bed"
#referencePath="/home/oabousaada/HPhasing/WAxNA/DeNovo/Reference/WAxNA.deNovo.fasta"
#minQ=1 #0.01
#minMQ=0 #40
#minAln=0.5
#outputPath="/home/oabousaada/HPhasing/WAxNA/VariantCalls/MinION/WAxNA.hetPositions.MinIONxSNPs.tsv"

#Initializing QScores

#P = 10-Q/10

QScores={"!":0,'"':1,"#":2,"$":3,"%":4,"&":5,"'":6,"(":7,")":8,"*":9,"+":10,",":11,"-":12,".":13,"/":14,"0":15,"1":16,"2":17,"3":18,"4":19,"5":20,"6":21,"7":22,"8":23,"9":24,":":25,";":26,"<":27,"=":28,">":29,"?":30,"@":31,"A":32,"B":33,"C":34,"D":35,"E":36,"F":37,"G":38,"H":39,"I":40,"J":41,"K":42} #Looks like the score is never 0 and that it's never over 30

PScores={}

for ASCII, Q in QScores.items():
    PScores[ASCII]=10**(-Q/10)

#Here it takes the different bedlines and it tries to apply the 0-base 1-base conversion maybe? #We'll only print things in these intervals
variantSet=set()

bedFile=open(bedPath,"r")

for line in bedFile:
    line=line.strip("\n").split("\t")
    contig=line[0]
    start=str(int(line[1])+1) #The bedFile is in 0-base, we want to stay in 1-base
    variantSet.add(contig+":"+start)

bedFile.close()

#This takes all the samfile lines
samLines=[]

samFile=open(samPath,"r")
for line in samFile:
    line=line.strip("\n").split("\t")
    samLines.append(line)

samFile.close()

readDict={}

for line in samLines:
    if line[0] not in readDict:
        readDict[line[0]]=[]
    else:
        i=1
        b=False
        while b==False:
            if line[0]+"_"+str(i) not in readDict:
                readDict[line[0]+"_"+str(i)]=[]
                line[0]=line[0]+"_"+str(i)
                b=True
            i+=1

#For some stats
samSkips=0
allAlnPcts=[]

i=0

for line in samLines: #For each sam record
    #Initializing a lot of things
    readName=line[0]
    ctgName=line[2]
    refPosition=int(line[3]) #SamFiles are aleady in 1-base
    MQ=int(line[4])
    cigar=line[5]
    sequence=line[9]
    QLine=line[10]
    readPosition=0 #The base of readPosition is irrelevant, it starts at 0 because python lists are 0-based (sequence[0])
    lineBases=[]
    matchStart=False
    averageScore=0
    for Q in QLine:
        averageScore+=QScores[Q]
    averageScore=int(averageScore/len(QLine))
    #
    #Okay so now we evaluate the cigar string and just slide down the SEQ and the reference simultaneously
    for match in re.finditer(cigarPartRE, cigar):
        n=int(match.group(1))
        alignmentType=match.group(2)
        if alignmentType=="S":
            readPosition+=n
        elif alignmentType in ["M", "=", "X"]:
            matchStart=True
            matches=[] #So, we know there might be some mismatches in here
            for x in range(n):
                matches.append((refPosition,sequence[readPosition]))
                refPosition+=1
                readPosition+=1
            for match in matches:
                pos=match[0]
                base=match[1]
                posName=ctgName+":"+str(pos)
                if posName in variantSet:
                    lineBases.append(posName+"="+base) #+"_"+str(QScores[Q])+"-"+str(averageScore)) #For debugging/info about which parameters can predict good reads
        elif alignmentType=="I":
            insertionSequence=""
            for x in range(n):
                insertionSequence+=sequence[readPosition]
                readPosition+=1
            pos=refPosition
            posName=ctgName+":"+str(pos)
            if posName in variantSet:
                #lineBases.append(posName+"="+insertionSequence) #+"_"+str(insertionScore)+"-"+str(averageScore)) #For debugging/info about which parameters can predict good reads
                pass #WE IGNORE THESE FOR NOW
        elif alignmentType=="D":
            if matchStart:
                for x in range(n):
                    matches.append((refPosition,"D",averageScore,averageScore))
                    refPosition+=1
                for match in matches:
                    pos=match[0]
                    base=match[1] #This will automatically be "D" for deletion
                    posName=ctgName+":"+str(pos)
                    if posName in variantSet:
                        #lineBases.append(posName+"="+base) #+"_"+str(Q)+"-"+str(averageScore)) #For debugging/info about which parameters can predict good reads
                        pass
    readDict[readName]=lineBases
    i+=1
    if i%(int(len(samLines)/100))==0:
        print(int(i/len(samLines)*100)+1,"%")

#skippedPct=(samSkips/len(samLines))*100

#print("There were",samSkips,"reads skipped because their alignments were smaller than the",minAln*100,"% minimum.\nThat's",round(skippedPct,2),"% of alignments in the sam file.")
#You can do something with allAlnPcts in order to determine if you want to change that percentage for whatever reason.

outputText=""

for readName, values in readDict.items():
    lineContents=readName+"\t"+"\t".join(values)+"\n"
    outputText+=lineContents

outputFile=open(outputPath,"w")
outputFile.write(outputText)
outputFile.close()

import os
import subprocess
import re
import gzip
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from plotnine import *

def longReadMapping(strainName,longReads,reference,outputFolder,flag,longReadPlatform,threads):
    outputLog=""
    #Map with ngmlr
    samFile=os.path.join(outputFolder,strainName+".sam")
    p=subprocess.run(["ngmlr","-t",threads,"-x",longReadPlatform,"-r",reference,"-q",longReads,"-o",samFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["ngmlr","-t",threads,"-x",longReadPlatform,"-r",reference,"-q",longReads,"-o",samFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Remove unmapped reads
    passSamFile=os.path.join(outputFolder,strainName+".pass.sam")
    p=subprocess.run(["samtools","view","-h","-t",reference,"-F",flag,samFile,"-o",passSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view","-t",reference,"-F",flag,samFile,"-o",passSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Sort sam file
    sortedHeaderSamFile=os.path.join(outputFolder,strainName+".sorted.header.sam")
    p=subprocess.run(["samtools","sort",passSamFile,"-o",sortedHeaderSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","sort",passSamFile,"-o",sortedHeaderSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #This looks stupid but is the safest way I could think of to get rid of the SAM header (it gets in the way of another function later)
    sortedSamFile=os.path.join(outputFolder,strainName+".sorted.sam")
    p=subprocess.run(["samtools","view",sortedHeaderSamFile,"-o",sortedSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view",sortedHeaderSamFile,"-o",sortedSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Cleaning up
    os.remove(samFile)
    os.remove(passSamFile)
    os.remove(sortedHeaderSamFile)

    return outputLog,"Long reads successfully mapped to reference"

def shortReadMapping(strainName,R1,R2,reference,outputFolder):
    outputLog=""
    #Mapping
    samFile=os.path.join(outputFolder, strainName+".sam")
    RGline="@RG\\tID:ID_"+strainName+"\\tLB:LB_"+strainName+"\\tPL:ILLUMINA\\tPU:PU_"+strainName+"\\tSM:SM_"+strainName
    p=subprocess.run(["bwa","mem","-M","-R", RGline, reference, R1, R2, "-o", samFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["bwa","mem","-M","-R", RGline, reference, R1, R2, "-o", samFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Organizing mapped reads
    bamFile=os.path.join(outputFolder,strainName+".bam")
    sortedBamFile=os.path.join(outputFolder,strainName+".sorted.bam")

    print("Removing quality 0 (multimapped) reads, turning to bam and sorting it")

    p=subprocess.run(["samtools","view","-bT",reference,"-q","1",samFile,"-o",bamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view","-bT",reference,"-q","1",samFile,"-o",bamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
    
    p=subprocess.run(["samtools", "sort", bamFile, "-o", sortedBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "sort", bamFile, "-o", sortedBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "index", sortedBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "index", sortedBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #GATK cleaning
    MDsortedBamFile=os.path.join(outputFolder,strainName+".sorted.MD.bam")
    p=subprocess.run(["gatk","MarkDuplicates", "--REMOVE_DUPLICATES", "true", "-I",sortedBamFile,"-O",MDsortedBamFile,"-M", "/dev/null"],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["gatk","MarkDuplicates", "--REMOVE_DUPLICATES", "true", "-I",sortedBamFile,"-O",MDsortedBamFile,"-M", "/dev/null"])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Finalizing
    finalBamFile=os.path.join(outputFolder,strainName+".final.bam")
    p=subprocess.run(["samtools", "sort", MDsortedBamFile, "-o", finalBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "sort", MDsortedBamFile, "-o", finalBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "index", finalBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "index", finalBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    try:
        #Cleaning up
        os.remove(samFile)
        os.remove(bamFile)
        os.remove(sortedBamFile)
        os.remove(MDsortedBamFile)
    except:
        return outputLog,"ERROR: Short reads unsuccessfully mapped to reference"

    return outputLog,"Short reads successfully mapped to reference"

def getShortReadPositions(hetVCFName,shortReadPosFileName):
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

    #Writing the short read het positions to a file
    shortReadSNPText=""

    for position in hetVCFPositions:
        shortReadSNPText+=position+"\n"

    shortReadSNPFile=open(shortReadPosFileName,"w")

    shortReadSNPFile.write(shortReadSNPText)

    shortReadSNPFile.close()

    return "Determined positions of heterozygous SNPs based on short read data"

def assignLongReadToSNPs(samPath,bedPath,referencePath,minQ,minMQ,minAln,outputPath):
    ###########
    #Resources#
    ###########

    #Cigar string:
    #http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
    #http://zenfractal.com/2013/06/19/playing-with-matches/
    #https://davetang.org/wiki/tiki-index.php?page=SAM

    #QScores
    #https://medium.com/@robertopreste/phred-quality-score-2837415f0af

    #To consider:
    #Currently, softclips, insertions and padding worsen the alignment score, you might want to change that

    cigarPartRE='(\d+)([DHIMNPSX=])'

    #Initializing QScores

    #P = 10-Q/10

    QScores={}
    for i in range(0,200):
        QScores[chr(33+i)]=i

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
        refPosition=int(line[3]) #SamFiles are already in 1-base
        MQ=int(line[4])
        cigar=line[5]
        sequence=line[9]
        QLine=line[10]
        readPosition=0 #The base of readPosition is irrelevant, it starts at 0 because python lists are 0-based (sequence[0])
        lineBases=[]
        matchStart=False
        averageScore=0
        for Q in QLine:
            try:
                averageScore+=QScores[Q]
            except:
                QScores[Q]=ord(Q)-33
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

    return "Reduced long reads to heterozygous SNPs based on short read information"

def longReadValidation(longReadFilePath,minCov,minRatio,minTrioCov,validatedSNPAssignmentsFile,contextDepthsFile):
    longReadFile=open(longReadFilePath,"r")

    longReads=[]

    for read in longReadFile:
        longReads.append(read.strip("\n").split("\t"))

    longReadFile.close()

    contextDepths={}
    distance=2

    for line in longReads:
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

    for line in longReads:
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


    longReadSets={}
    longReadTuples={}


    ###########################################################################################
    #You can't just do this step before you remove reads based on coverage, that's ridiculous.#
    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
    ###########################################################################################

    #Only keeping long reads that have unique combinations of SNPs 
    uniqueLongReadSets=set()
    uniqueLongReadIDs=set()

    curLen=len(uniqueLongReadSets)

    for read in longReads:
        readSet=frozenset(read[1:])
        altruisticSet=[]
        for SNP in readSet:
            altruisticSet.append(SNP)
        readSet=frozenset(altruisticSet)
        uniqueLongReadSets.add(readSet)
        if len(uniqueLongReadSets)>curLen:
            curLen+=1
            longReadSets[read[0]]=readSet
            longReadTuples[read[0]]=tuple(readSet)
            uniqueLongReadIDs.add(read[0])

    ###########################################################################################
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
    #You can't just do this step before you remove reads based on coverage, that's ridiculous.#
    ###########################################################################################

    SNPAssignments=longReadSets

    #Printing it all to a file
    SNPAssignmentsText=""

    for longRead, assignedSNPs in SNPAssignments.items():
        if assignedSNPs!=set({''}):
            SNPAssignmentsText+=longRead+"\t"+"\t".join(assignedSNPs)+"\n"

    SNPAssignmentFile=open(validatedSNPAssignmentsFile,"w")
    SNPAssignmentFile.write(SNPAssignmentsText)
    SNPAssignmentFile.close()

    return "Pre-processed long reads and generated context depth file"

def generateLongReadFastQFiles(haplotigReadNameFilePath,longReadFastQFilePath,outputPath):
    haplotigReadDict={}
    haplotigReadNameFile=open(haplotigReadNameFilePath,"r")
    for line in haplotigReadNameFile:
        line=line.strip("\n").split("\t")
        if line[0] not in haplotigReadDict:
            haplotigReadDict[line[0]]=[line[1]]
        else:
            haplotigReadDict[line[0]].append(line[1])

    haplotigReadNameFile.close()

    gzipBool=True
    
    longReadData={}

    with gzip.open(longReadFastQFilePath, 'r') as fh:
        try:
            longReadFastQFile=gzip.open(longReadFastQFilePath,"rt")
        except gzip.BadGzipFile:
            longReadFastQFile=open(longReadFastQFilePath,"r")
    i=0
    for line in longReadFastQFile:
        line=line.strip("\n")
        if i%4==0:
            readName=line
            longReadData[readName]=[]
        else:
            longReadData[readName].append(line)
        i+=1
    longReadFastQFile.close()

    for haplotigName, reads in haplotigReadDict.items():
        haplotigFastQText=""
        for read in reads:
            read="@"+read
            readText=read+"\n"+"\n".join(longReadData[read])+"\n"
            haplotigFastQText+=readText
        haplotigFile=gzip.open(outputPath+haplotigName+".fastq.gz",'wt')
        haplotigFile.write(haplotigFastQText)
        haplotigFile.close()

    return "Successfully generated phased FastQ files"

def loadFile(filePath):
	openFile=open(filePath,"r")
	fileContents=[]
	for line in openFile:
		line=line.strip("\n").split("\t")
		fileContents.append(line)
	openFile.close()
	return fileContents

def simplifyDataVis(dataVisPath,simpleOutPath,distance):
    fileContents=loadFile(dataVisPath)
    fileDict={}

    for line in fileContents:
        if line[0] not in fileDict:
            fileDict[line[0]]=[line]
        else:
            fileDict[line[0]].append(line)

    newLines=[]

    for haplotig, SNPLines in fileDict.items():
        startPos=int(SNPLines[0][1])
        endPos=int(SNPLines[0][1])+1
        for line in SNPLines:
            if int(line[1])>endPos+distance:
                newLines.append([haplotig,str(startPos),str(endPos),line[2],line[3]])
                startPos=int(line[1])
                endPos=int(line[1])+1
            else:
                endPos=int(line[1])
        newLines.append([haplotig, str(startPos), str(endPos), line[2], line[3]])

    outFileText=""
    for line in newLines:
        outFileText+="\t".join(line)+"\n"

    outFile=open(simpleOutPath,"w")
    outFile.write(outFileText)
    outFile.close()

def generateDataVis(dataVisPath,outPath):
    outputPNG=outPath+"dataVis.png"
    outputPDF=outPath+"dataVis.pdf"
    outputSVG=outPath+"dataVis.svg"

    tbl = pd.read_csv(dataVisPath, sep="\t", header=None)
    tbl.columns = ["contigName", "startPos", "endPos", "chr", "yValue"]

    #g=(ggplot(tbl,aes(y=tbl["yValue"],yend=tbl["yValue"],x=tbl["startPos"],xend=tbl["endPos"],color='factor(yValue)'))+geom_segment(size=1.5)+theme(legend_position="none")+theme(panel_grid_minor=element_blank())+facet_wrap("~chr",scales = "free")+theme(axis_title_y=element_blank(),axis_text_y=element_blank(),axis_ticks_major_y=element_blank())+xlab("Distance from start of genome (bp)")+theme(subplots_adjust={'hspace': 0.7}))
    
    g=(ggplot(tbl,aes(y=tbl["yValue"],yend=tbl["yValue"],x=tbl["startPos"],xend=tbl["endPos"],color='factor(yValue)'))+geom_segment(size=1.5)+theme(legend_position="none")+theme(panel_grid_minor=element_blank())+theme(axis_title_y=element_blank(),axis_text_y=element_blank(),axis_ticks_major_y=element_blank())+xlab("Distance from start of genome (bp)")+theme(subplots_adjust={'hspace': 0.7}))

    ggsave(g,filename=outputSVG,width=18, height=10)
    ggsave(g,filename=outputPNG,width=18, height=10)
    ggsave(g,filename=outputPDF,width=18, height=10)

    return "Generated plots"


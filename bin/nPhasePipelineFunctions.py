import os
import subprocess
import re
import gzip
import sortedcontainers
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from plotnine import *
import multiprocessing

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
    p=subprocess.run(["samtools","view","-h","-t",reference,"-@",threads,"-F",flag,samFile,"-o",passSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view","-t",reference,"-@",threads,"-F",flag,samFile,"-o",passSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Sort sam file
    sortedHeaderSamFile=os.path.join(outputFolder,strainName+".sorted.header.sam")
    p=subprocess.run(["samtools","sort",passSamFile,"-@",threads,"-o",sortedHeaderSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","sort",passSamFile,"-@",threads,"-o",sortedHeaderSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #This looks stupid but is the safest way I could think of to get rid of the SAM header (it gets in the way of another function later)
    sortedSamFile=os.path.join(outputFolder,strainName+".sorted.sam")
    p=subprocess.run(["samtools","view",sortedHeaderSamFile,"-@",threads,"-o",sortedSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view",sortedHeaderSamFile,"-@",threads,"-o",sortedSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    try:
        #Cleaning up
        os.remove(samFile)
        os.remove(passSamFile)
        os.remove(sortedHeaderSamFile)
    except:
        return outputLog,"Long reads not mapped to reference"

    return outputLog,"Long reads successfully mapped to reference"

def shortReadMapping(strainName,R1,R2,reference,outputFolder,threads):
    outputLog=""

    #Mapping
    samFileNoRGPath=os.path.join(outputFolder, strainName+".noRG.sam")
    samFileNoRG=open(samFileNoRGPath,"w")
    p=subprocess.run(["bwa","mem","-t",threads,"-M", reference, R1, R2],stderr=subprocess.PIPE,stdout=samFileNoRG, universal_newlines=True)
    samFileNoRG.close()
    outputLog+="COMMAND: "+" ".join(["bwa","mem","-t",threads,"-M", reference, R1, R2, ">", samFileNoRGPath])+"\n\n"
    if p.stderr!="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT is in "+samFileNoRGPath+"\n\n"

    #Adding read groups
    samFile=os.path.join(outputFolder, strainName+".sam")
    p=subprocess.run(["gatk", "AddOrReplaceReadGroups","-I",samFileNoRGPath,"-O",samFile,"-RGID","ID_"+strainName, "-RGLB","LB_"+strainName, "-RGPL","ILLUMINA", "-RGPU","PU_"+strainName, "-RGSM","SM_"+strainName],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["gatk", "AddOrReplaceReadGroups","-I",samFileNoRGPath,"-O",samFile,"-RGID","ID_"+strainName, "-RGLB","LB_"+strainName, "-RGPL","ILLUMINA", "-RGPU","PU_"+strainName, "-RGSM","SM_"+strainName])+"\n\n"
    if p.stderr!="" or p.stdout!="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Organizing mapped reads
    bamFile=os.path.join(outputFolder,strainName+".bam")
    sortedBamFile=os.path.join(outputFolder,strainName+".sorted.bam")

    print("Removing quality 0 (multimapped) reads, turning to bam and sorting it")

    p=subprocess.run(["samtools","view","-@",threads,"-bT",reference,"-q","1",samFile,"-o",bamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view","-@",threads,"-bT",reference,"-q","1",samFile,"-o",bamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "sort", bamFile,"-@",threads, "-o", sortedBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "sort", bamFile,"-@",threads, "-o", sortedBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "index","-@",threads, sortedBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "index","-@",threads, sortedBamFile])+"\n\n"
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
    p=subprocess.run(["samtools", "sort", MDsortedBamFile, "-@",threads,"-o", finalBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "sort", MDsortedBamFile, "-@",threads,"-o", finalBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "index","-@",threads, finalBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "index","-@",threads, finalBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    try:
        #Cleaning up
        os.remove(samFileNoRGPath)
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
            chr=line[0].replace(":","_")
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
                    #   hetVCFPositions.add(chr+":"+str(curPos))
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
        contig=line[0].replace(":","_")
        start=str(int(line[1])+1) #The bedFile is in 0-base, we want to stay in 1-base
        variantSet.add(contig+":"+start)

    bedFile.close()

    #This takes all the samfile lines
    samLines=[]

    samFile=open(samPath,"r")
    for line in samFile:
        line=line.strip("\n").split("\t")
        if line[0][0]!="@":
            line[2]=line[2].replace(":","_")
            samLines.append(line)

    samFile.close()

    manager = multiprocessing.Manager()

    readDict = manager.dict()

    with multiprocessing.Pool(processes=8) as pool: #Consider batching lines
        pool.map(assignLongReadToSNPsParallelHelper, [[line,readDict,QScores,variantSet] for line in samLines])

    outputText=""

    for readName, values in readDict.items():
        lineContents=readName+"\t"+"\t".join(values)+"\n"
        outputText+=lineContents

    outputFile=open(outputPath,"w")
    outputFile.write(outputText)
    outputFile.close()

    return "Reduced long reads to heterozygous SNPs based on short read information"

def assignLongReadToSNPsParallelHelper(variables):
    line=variables[0]
    rD=variables[1]
    QScores=variables[2]
    variantSet=variables[3]
    assignLongReadToSNPsParallel(line,rD,QScores,variantSet)

def assignLongReadToSNPsParallel(line,readDict,QScores,variantSet):
    cigarPartRE = '(\d+)([DHIMNPSX=])'
    if line[0] not in readDict:
        readDict[line[0]]=[]
    else:
        i=1
        b=False
        while b==False:
            if line[0]+"_"+str(i) not in readDict:
                readDict[line[0]+"_VCSTTP_"+str(i)]=[] #Very Clean Solution To This Problem
                line[0]=line[0]+"_VCSTTP_"+str(i)
                b=True
            i+=1

    #Initializing a lot of things
    sequence=line[9]
    if sequence=="*": #Only seen this happen when the MAPQ is 0 for a secondary alignment
        return
    readName=line[0]
    ctgName=line[2]
    refPosition=int(line[3]) #SamFiles are already in 1-base
    MQ=int(line[4])
    cigar=line[5]
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
                        pass #WE IGNORE THESE FOR NOW
    readDict[readName]=lineBases

def longReadValidation(longReadFilePath,minCov,minRatio,minTrioCov,validatedSNPAssignmentsFile,contextDepthsFile):
    longReadFile=open(longReadFilePath,"r")

    longReads=[]

    for read in longReadFile:
        longReads.append(read.strip("\n").split("\t"))

    longReadFile.close()

    contextDepths={}
    batchSize=3000
    batches=[]
    currentBatch=[]
    i=0
    for longRead in longReads:
        currentBatch.append(longRead)
        i+=1
        if i>batchSize:
            batches.append(currentBatch)
            currentBatch=[]
            i=0

    if currentBatch!=[]:
        batches.append(currentBatch)

    with multiprocessing.Pool(processes=8) as pool:
        contextDepthDicts=pool.map(updateContextDepths, [batchData for batchData in batches])
        for contextDepthDict in contextDepthDicts:
            for SNP, SNPContexts in contextDepthDict.items():
                contextDepths.setdefault(SNP,{})
                for SNPContext, SNPContextCount in SNPContexts.items():
                    contextDepths[SNP].setdefault(SNPContext,0)
                    contextDepths[SNP][SNPContext]+=SNPContextCount

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

def updateContextDepths(longReadBatch):
    localContextDepths={}
    distance=2
    for longRead in longReadBatch:
        SNPs=longRead[1:]
        for SNP in SNPs:
            i=SNPs.index(SNP)
            itemList=[]
            for j in range(i-distance,i+distance+1):
                if j>=0 and j<len(SNPs):
                    itemList.append(SNPs[j])
            if SNP!='':
                frozenList=frozenset(itemList)
                localContextDepths.setdefault(SNP,{})
                localContextDepths[SNP].setdefault(frozenList,0)
                localContextDepths[SNP][frozenList]+=1
    return localContextDepths

def getAvgFastQScore(fastQScoreStr):
    fastQScoreStr=fastQScoreStr.strip()
    totalScore=0
    QScores={}
    for baseQ in fastQScoreStr:
        if baseQ not in QScores:
            QScores[baseQ]=ord(baseQ)-33
        totalScore+=QScores[baseQ]
    avgScore=totalScore/len(fastQScoreStr)
    return avgScore

def generateLongReadFastQFiles(haplotigReadNameFilePath,longReadFastQFilePath,outputPath): #Parallelize in a way that can save memory
    ##clusterStatsText=""
    haplotigReadDict={}
    allAllowedReads=set()
    haplotigReadNameFile=open(haplotigReadNameFilePath,"r")
    for line in haplotigReadNameFile:
        line=line.strip("\n").split("\t")
        haplotigName=line[0]
        readName="@"+line[1]
        readName=readName.split("_VCSTTP_")[0]
        haplotigReadDict.setdefault(readName,set())
        haplotigReadDict[readName].add(haplotigName)
        allAllowedReads.add(readName)

    haplotigReadNameFile.close()


    ##longReadMetaData={}

    try:
        longReadFastQFile=gzip.open(longReadFastQFilePath,"rt")
        longReadFastQFile.readline()
        longReadFastQFile.close()
        longReadFastQFile=gzip.open(longReadFastQFilePath, "rt")
    except gzip.BadGzipFile:
        longReadFastQFile=open(longReadFastQFilePath,"r")

    fastQReadQueue = multiprocessing.Queue()

    fastQWriter_p = multiprocessing.Process(target=fastQWriter, args=((outputPath),(fastQReadQueue),))
    fastQWriter_p.daemon = True
    fastQWriter_p.start()

    i=0
    j=0
    putCounter=0
    readData=[]
    readDataBatch={}
    batchSize=2500
    for line in longReadFastQFile:
        line=line.strip("\n")
        if i%4==0:
            if readData!=[]:
                for haplotigName in haplotigNames:
                    readDataBatch.setdefault(haplotigName,[])
                    readDataBatch[haplotigName].append(readData)
                    j+=1
                if j>batchSize:
                    j=0
                    fastQReadQueue.put(readDataBatch)
                    putCounter+=1
                    print(putCounter*batchSize)
                    readDataBatch={}
            line=line.split()
            readName=line[0]
            haplotigNames=[]
            readData=[]
            if readName in allAllowedReads:
                for haplotigName in haplotigReadDict[readName]:
                    haplotigNames.append(haplotigName)
                readData=[readName]
        else:
            if readName in allAllowedReads:
                ##if i%4==3:
                    ##longReadMetaData[readName]=getAvgFastQScore(line)
                readData.append(line)
        i+=1
    longReadFastQFile.close()

    print(putCounter * batchSize + j)

    if readDataBatch!={}:
        fastQReadQueue.put(readDataBatch)

    fastQReadQueue.put("end")

    fastQWriter_p.join()

    ##clusterStatsFile=open(outputPath+"clusterStats.tsv", 'w')
    ##clusterStatsFile.write(clusterStatsText)
    ##clusterStatsFile.close()

    return "Successfully generated phased FastQ files"

def fastQWriter(outputPath,queue):
    while True:
        readDataBatch=queue.get()
        if readDataBatch == "end":
            return
        for haplotigName, haplotigReads in readDataBatch.items():
            haplotigFastQText=""
            for readData in haplotigReads:
                readName=readData[0]
                readString=readName+"\n"+"\n".join(readData[1:])+"\n"
                haplotigFastQText+=readString
            haplotigFile=gzip.open(outputPath+haplotigName+".fastq.gz",'at')
            haplotigFile.write(haplotigFastQText)
            haplotigFile.close()

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

def generatePhasingVis(dataVisPath,outPath):
    #Figure out how supported (as a %) each base is for each position in the cluster

    outputPNG=outPath+"phasedVis.png"
    outputPDF=outPath+"phasedVis.pdf"
    outputSVG=outPath+"phasedVis.svg"

    tbl=pd.read_csv(dataVisPath,sep="\t",header=None)
    tbl.columns=["contigName","startPos","endPos","chr","yValue"]

    g=(ggplot(tbl,aes(y="contigName",yend="contigName",x="startPos",xend="endPos",color='contigName'))+geom_segment(size=1.5)+theme(legend_position="none")+theme(panel_grid_minor=element_blank())+facet_wrap("~chr",scales="free")+theme(axis_title_y=element_blank(),axis_text_y=element_blank(),axis_ticks_major_y=element_blank())+xlab("Position (bp)")+theme(subplots_adjust={'hspace':0.7}))

    ggsave(g,filename=outputSVG,width=18,height=10)
    ggsave(g,filename=outputPNG,width=18,height=10)
    ggsave(g,filename=outputPDF,width=18,height=10)

    return "Generated phased plots"

def generateCoverageVis(dataVisPath,outPath):
    #Figure out how covered (as X) each heterozygous position is for each cluster

    outputPNG=outPath+"coverageVis.png"
    outputPDF=outPath+"coverageVis.pdf"
    outputSVG=outPath+"coverageVis.svg"

    tbl=pd.read_csv(dataVisPath,sep="\t",header=None,na_values="NA")
    tbl.columns=["haplotigName","chr","pos","coverage"]

    g=(ggplot(tbl,aes(x="pos",y="coverage",color='haplotigName'))+geom_line(size=1.5)+theme(legend_position="none")+theme(panel_grid_minor=element_blank())+facet_wrap("~chr",scales="free")+ylab("Coverage (X)")+xlab("Position (bp)")+theme(subplots_adjust={'hspace': 0.7}))

    ggsave(g,filename=outputSVG,width=18,height=10)
    ggsave(g,filename=outputPNG,width=18,height=10)
    ggsave(g,filename=outputPDF,width=18,height=10)

    return "Generated coverage plots"

def generateDiscordanceVis(dataVisPath,outPath): #This function is very slow with plotnine for some reason. Should replace this in the long run.
    #Figure out how supported (as a %) each base is for each position in the cluster

    #Takes several minutes to output each graph and violin plots have an unreasonable memory constraint using plotnine.
    #TODO: Fix this
    #TEMP: Moved this to very last part of the pipeline so that everything else can be generated smoothly first.
    #TEMP: Commented out two of the less useful file formats to save time.

    #outputPNG=outPath+"discordanceVis.png"
    #outputPDF=outPath+"discordanceVis.pdf"
    outputSVG=outPath+"discordanceVis.svg"

    tbl=pd.read_csv(dataVisPath,sep="\t",header=None,comment="#")
    tbl.columns=["cluster","chr","position","base","frequency","coverage"]

    g=(ggplot(tbl,aes(y="frequency",x='cluster'))+facet_wrap("~chr",scales="free")+geom_violin(scale="width")+coord_flip()+ylim(0,1)+xlab("Cluster name")+ylab("Allele frequency (%)")+ggtitle("Allele frequency by cluster (all clusters have same max width)")+theme(subplots_adjust={'wspace': 0.25}))

    #ggsave(g,filename=outputPNG,width=18,height=10)
    #ggsave(g,filename=outputPDF,width=18,height=10)
    ggsave(g,filename=outputSVG,width=18,height=10)

    return "Generated discordance plots"

def generateDiscordance(clusterReadFilePath,readFilePath,outPath):
    #Figure out how supported (as a %) each base is for each position in the cluster

    #Identify warnings/sort by how bad & output a text file

    #Output a plot
    clusterDict={} #All reads associated with clusters

    clusterFile=open(clusterReadFilePath,"r")

    for line in clusterFile:
        line=line.strip("\n").split("\t")
        if line[0] not in clusterDict:
            clusterDict[line[0]]=set()
        clusterDict[line[0]].add(line[1])

    clusterFile.close()

    readDict={} #All SNPs associated with reads

    readFile=open(readFilePath,"r")

    for line in readFile:
        line=line.strip("\n").split("\t")
        if line[0] not in readDict:
            readDict[line[0]]=set()
        readDict[line[0]].update(line[1:])

    readFile.close()

    rawDict={} #All SNPs within clusters based on the reads involved

    for cluster,reads in clusterDict.items():
        rawDict[cluster]={}
        for read in reads:
            for SNP in readDict[read]:
                chr=SNP.split(":")[0]
                position=SNP.split(":")[1].split("=")[0]
                base=SNP.split("=")[1]
                if chr not in rawDict[cluster]:
                    rawDict[cluster][chr]={}
                if position not in rawDict[cluster][chr]:
                    rawDict[cluster][chr][position]={}
                if base not in rawDict[cluster][chr][position]:
                    rawDict[cluster][chr][position][base]=0
                rawDict[cluster][chr][position][base]+=1 #Count the number for each base

    allFrequencies=[]

    for cluster,chrs in rawDict.items():
        for chr,positions in chrs.items():
            for position, bases in positions.items():
                baseCovs=bases.values()
                totalCov=sum(baseCovs)
                for base, coverage in bases.items():
                    allFrequencies.append((cluster,chr,position,base,str(coverage/totalCov),str(totalCov)))

    cleanFullText="#cluster\tchr\tposition\tbase\tfrequency\tcoverage\n"

    for freq in allFrequencies:
        cleanFullText+="\t".join(freq)+"\n"

    cleanFullTextOutputFile=open(outPath,"w")
    cleanFullTextOutputFile.write(cleanFullText)
    cleanFullTextOutputFile.close()

    pass

def generateCoverage(clusterFilePath,readFilePath,outPath,windowSize):
    #Figure out how covered each read position is in each cluster

    totalCoverage=0
    uniquePositions=set()

    clusterDict={}

    clusterFile=open(clusterFilePath,"r")

    for line in clusterFile:
        line=line.strip("\n").split("\t")
        if line[0] not in clusterDict:
            clusterDict[line[0]]=set()
        clusterDict[line[0]].add(line[1])

    clusterFile.close()

    readDict={}

    readFile=open(readFilePath,"r")

    for line in readFile:
        line=line.strip("\n").split("\t")
        if line[0] not in readDict:
            readDict[line[0]]=set()
        readDict[line[0]].update(line[1:])

    readFile.close()

    rawDict={}

    for cluster,reads in clusterDict.items():
        rawDict[cluster]={}
        for read in reads:
            for SNP in readDict[read]:
                uniquePositions.add(SNP.split("=")[0])
                chr=SNP.split(":")[0]
                position=SNP.split(":")[1].split("=")[0]
                if chr not in rawDict[cluster]:
                    rawDict[cluster][chr]={}
                if position not in rawDict[cluster][chr]:
                    rawDict[cluster][chr][position]=0
                rawDict[cluster][chr][position]+=1
                totalCoverage+=1

    cleanFullText=""

    betterClusters={}

    for cluster,chrs in rawDict.items():
        betterClusters[cluster]={}
        for chr,positions in rawDict[cluster].items():
            betterClusters[cluster][chr]=[]
            for pos in positions:
                cov=str(rawDict[cluster][chr][pos])
                if chr not in betterClusters[cluster].keys():
                    betterClusters[cluster][chr]=[]
                betterClusters[cluster][chr].append((int(pos),int(cov)))
            sortedPositions=betterClusters[cluster][chr]
            sortedPositions.sort()
            betterClusters[cluster][chr]=sortedPositions

    windowDict={}

    for cluster,chrs in betterClusters.items():
        windowDict[cluster]={}
        for chr in chrs.keys():
            if chr not in windowDict[cluster]:
                windowDict[cluster][chr]={}
            sortedPositions=betterClusters[cluster][chr]
            sortedPositions.sort()
            minPos=0
            maxPos=int(sortedPositions[-1][0])+windowSize
            for i in range(minPos,maxPos,windowSize):
                window=[]
                for position in sortedPositions:
                    if int(position[0])>=i and int(position[0])<i+windowSize:
                        window.append(int(position[1]))
                if len(window)==0:
                    mean="NA"
                else:
                    mean=sum(window)/len(window)
                    if mean==0:
                        mean="NA"
                windowDict[cluster][chr][i+(windowSize/2)]=str(mean)

    cleanFullText=""

    for cluster,chrs in windowDict.items():
        for chr,positions in chrs.items():
            for pos,cov in positions.items():
                cleanFullText+=cluster+"\t"+chr+"\t"+str(pos)+"\t"+str(cov)+"\n"

    cleanFullTextOutputFile=open(outPath,"w")
    cleanFullTextOutputFile.write(cleanFullText)
    cleanFullTextOutputFile.close()

    minCov=int(0.05*(totalCoverage/len(uniquePositions)))

    return minCov


def giveMeFullData(clusters):
    clusterText=""
    sortedClusterLines=sortedcontainers.SortedList()
    clusterLines=[]
    for clusterName, cluster in clusters.items():
        for SNP in cluster:
            contig=SNP.split(":")[0]
            position=int(SNP.split(":")[1].split("=")[0])
            sortedClusterLines.add([position,clusterName,contig])
    i=1
    previouslySeen={sortedClusterLines[0][1]:i}
    for line in sortedClusterLines:
        if line[1] not in previouslySeen:
            i+=1
            previouslySeen[line[1]]=i
        clusterLines.append([line[1],line[0],line[2],previouslySeen[line[1]]])
    for line in clusterLines:
        clusterText+="\t".join([str(x) for x in line])+"\n"
    return clusterText

def filterCoverage(clusterReadFilePath,readFilePath,outPath,minCov):
    #Filter out how based on coverage

    clusterDict={} #All reads associated with clusters

    clusterFile=open(clusterReadFilePath,"r")

    for line in clusterFile:
        line=line.strip("\n").split("\t")
        if line[0] not in clusterDict:
            clusterDict[line[0]]=set()
        clusterDict[line[0]].add(line[1])

    clusterFile.close()

    readDict={} #All SNPs associated with reads

    readFile=open(readFilePath,"r")

    for line in readFile:
        line=line.strip("\n").split("\t")
        if line[0] not in readDict:
            readDict[line[0]]=set()
        readDict[line[0]].update(line[1:])

    readFile.close()

    rawDict={} #All SNPs within clusters based on the reads involved

    for cluster,reads in clusterDict.items():
        rawDict[cluster]={}
        for read in reads:
            for SNP in readDict[read]:
                chr=SNP.split(":")[0]
                position=SNP.split(":")[1].split("=")[0]
                base=SNP.split("=")[1]
                if chr not in rawDict[cluster]:
                    rawDict[cluster][chr]={}
                if position not in rawDict[cluster][chr]:
                    rawDict[cluster][chr][position]=0
                rawDict[cluster][chr][position]+=1 #Count the number for each base

    allFrequencies=[]

    for cluster,chrs in rawDict.items():
        for chr,positions in chrs.items():
            for position, coverage in positions.items():
                if coverage>minCov:
                    allFrequencies.append((cluster,chr,position,str(coverage)))

    cleanFullText="#cluster\tchr\tcoverage\tposition\n"

    for freq in allFrequencies:
        cleanFullText+="\t".join(freq)+"\n"

    cleanFullTextOutputFile=open(outPath,"w")
    cleanFullTextOutputFile.write(cleanFullText)
    cleanFullTextOutputFile.close()

    pass






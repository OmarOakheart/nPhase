import argparse
import sys
import os
import subprocess
import bin.nPhasePipelineFunctions as nPhaseFunctions
import bin.nPhase as phaseTool
import bin.nPhaseCleaning as cleaningTool

def updateLog(logFilePath,logText):
    logFile=open(logFilePath,"a")
    logFile.write(logText)
    logFile.close()


def nPhasePipeline(args):
    #Create folder structure
    basePath=os.path.join(args.outputFolder,args.strainName)
    mappedShortReadPath=os.path.join(basePath,"Mapped","shortReads")
    variantCalledShortReadPath=os.path.join(basePath,"VariantCalls","shortReads")
    mappedLongReadPath=os.path.join(basePath,"Mapped","longReads")
    variantCalledLongReadPath=os.path.join(basePath,"VariantCalls","longReads")
    overlapPath=os.path.join(basePath,"Overlaps")
    phasedPath=os.path.join(basePath,"Phased")
    phasedFastqPath=os.path.join(basePath,"Phased","FastQ")
    datavisFolderPath=os.path.join(basePath,"Phased","Plots")
    logPath=os.path.join(basePath,"Logs") #Currently not producing any logs

    allPaths=[basePath,mappedShortReadPath,variantCalledShortReadPath,mappedLongReadPath,variantCalledLongReadPath,overlapPath,phasedPath,phasedFastqPath,datavisFolderPath,logPath]

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    readmePath=os.path.join(basePath,"Readme.txt")
    fullLogPath=os.path.join(basePath,"Logs","fullLog.txt")

    ########################
    #Pre-process reference #
    ########################

    #Make sure the reference is indexed (might be missing a process)
    if not os.path.isfile(args.reference+".fai"):
        p=subprocess.run(["samtools","faidx",args.reference],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["samtools","faidx",args.reference])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)

    if not os.path.isfile(args.reference+".bwt"):
        p=subprocess.run(["bwa","index",args.reference],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["bwa","index",args.reference])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)

    if not os.path.isfile(args.reference+".dict"):
        p=subprocess.run(["gatk","CreateSequenceDictionary","-R",args.reference],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["gatk","CreateSequenceDictionary","-R",args.reference])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)

    ########################
    #Pre-process long reads#
    ########################

    #Map long reads to reference
    splitReadFlag="260" #This flag allows split reads
    outputLog,systemMessage=nPhaseFunctions.longReadMapping(args.strainName,args.longReadFile,args.reference,mappedLongReadPath,splitReadFlag,args.longReadPlatform,args.threads)
    print(systemMessage)
    updateLog(fullLogPath,outputLog)


    #########################
    #Pre-process short reads#
    #########################

    #Map short reads to reference
    outputLog,systemMessage=nPhaseFunctions.shortReadMapping(args.strainName,args.shortReadFile_R1,args.shortReadFile_R2,args.reference,mappedShortReadPath,args.threads)
    print(systemMessage)
    updateLog(fullLogPath,outputLog)

    #Variant call short reads and select SNPs only
    estimatedPloidy="2" #This argument is required for GATK's variant calling and we can expect most SNPs to only have two alleles anyway

    shortReadBam=os.path.join(mappedShortReadPath,args.strainName+".final.bam")
    shortReadVCF=os.path.join(variantCalledShortReadPath,args.strainName+".vcf")
    p=subprocess.run(["gatk","HaplotypeCaller","-R",args.reference,"-ploidy",estimatedPloidy,"-I",shortReadBam,"-O",shortReadVCF],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    logText="COMMAND: "+" ".join(["gatk","HaplotypeCaller","-R",args.reference,"-ploidy",estimatedPloidy,"-I",shortReadBam,"-O",shortReadVCF])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
    updateLog(fullLogPath,logText)

    shortReadSNPsVCF=os.path.join(variantCalledShortReadPath,args.strainName+".SNPs.vcf")
    p=subprocess.run(["gatk","SelectVariants","-R",args.reference,"--variant",shortReadVCF,"-O",shortReadSNPsVCF,"--select-type-to-include","SNP"],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    logText="COMMAND: "+" ".join(["gatk","SelectVariants","-R",args.reference,"--variant",shortReadVCF,"-O",shortReadSNPsVCF,"--select-type-to-include","SNP"])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
    updateLog(fullLogPath,logText)


    #Extract heterozygous positions from VCF file
    #Ugh this is a whole thing.

    SNPVCFFilePath=os.path.join(variantCalledShortReadPath,args.strainName+".SNPs.vcf")
    hetSNPVCFFilePath=os.path.join(variantCalledShortReadPath,args.strainName+".hetSNPs.vcf")

    hetSNPVCFText=""

    SNPVCFFile=open(SNPVCFFilePath,"r")
    for line in SNPVCFFile:
        if "#" in line:
            hetSNPVCFText+=line
        elif "AF=1.00" not in line:
            line=line.strip("\n")
            line=line.replace(";","\t")
            line=line.split("\t")
            line=[line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7]]
            hetSNPVCFText+="\t".join(line)+"\n"

    SNPVCFFile.close()

    hetSNPVCFFile=open(hetSNPVCFFilePath,"w")
    hetSNPVCFFile.write(hetSNPVCFText)
    hetSNPVCFText="" #Just saving some memory since this isn't in a function
    hetSNPVCFFile.close()

    print("Identified heterozygous SNPs in short read VCF")

    shortReadPositionsOutputFilePath=os.path.join(mappedShortReadPath,args.strainName+".hetSNPs.positions.tsv")
    nPhaseFunctions.getShortReadPositions(hetSNPVCFFilePath,shortReadPositionsOutputFilePath)

    shortReadSNPsBedFilePath=os.path.join(mappedShortReadPath,args.strainName+".hetSNPs.bed")

    shortReadSNPsBedText=""

    shortReadPositionsOutputFile=open(shortReadPositionsOutputFilePath,"r")
    for line in shortReadPositionsOutputFile:
        if "#" in line:
            shortReadSNPsBedText+=line
        elif "AF=1.00" not in line:
            line=line.strip("\n")
            line=line.replace(":","\t")
            line=line.split("\t")
            line=[line[0],str(int(line[1])-1),line[1]]
            shortReadSNPsBedText+="\t".join(line)+"\n"

    shortReadPositionsOutputFile.close()

    shortReadSNPsBedFile=open(shortReadSNPsBedFilePath,"w")
    shortReadSNPsBedFile.write(shortReadSNPsBedText)
    shortReadSNPsBedText="" #Just saving some memory since this isn't in a function
    shortReadSNPsBedFile.close()

    print("Extracted heterozygous SNP info based on short read VCF")

    #Reduce long reads to their heterozygous SNPs
    cleanLongReadSamFile=os.path.join(mappedLongReadPath,args.strainName+".sorted.sam")
    minQ=0.01   # Currently
    minMQ=0     # Not
    minAln=0.5  # Used
    longReadPositionNTFile=os.path.join(variantCalledLongReadPath,args.strainName+".hetPositions.SNPxLongReads.tsv")
    nPhaseFunctions.assignLongReadToSNPs(cleanLongReadSamFile,shortReadSNPsBedFilePath,args.reference,minQ,minMQ,minAln,longReadPositionNTFile)

    #Only keep the longest reads and get rid of duplicate heterozygous SNP profiles, compensate by keeping context coverage information.
    minCov=0     # Currently
    minRatio=0   # Not
    minTrioCov=0 # Used
    validatedSNPAssignmentsFile=os.path.join(variantCalledLongReadPath,args.strainName+".hetPositions.SNPxLongReads.validated.tsv")
    contextDepthsFile=os.path.join(overlapPath,args.strainName+".contextDepths.tsv")
    nPhaseFunctions.longReadValidation(longReadPositionNTFile,minCov,minRatio,minTrioCov,validatedSNPAssignmentsFile,contextDepthsFile)

    #Run everything through nPhase
    phaseTool.nPhase(validatedSNPAssignmentsFile,args.strainName,contextDepthsFile,phasedPath,basePath,args.reference,args.minSim,args.minOvl,args.minLen,args.maxID)

    readmeText="\nPhased files can be found at "+phasedPath+"\nThe *_variants.tsv file contains information on the consensus heterozygous variants present in each predicted haplotig.\nThe *_clusterReadNames.tsv file contains information on the reads which comprise each cluster."
    print(readmeText)
    updateLog(readmePath,readmeText)

    #Simplify datavis
    dataVisPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_phasedDataFull.tsv")
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_phasedDataSimple.tsv")
    nPhaseFunctions.simplifyDataVis(dataVisPath,simpleOutPath,1000)

    #Generate plots
    datavisPath=os.path.join(datavisFolderPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_")
    nPhaseFunctions.generatePhasingVis(simpleOutPath,datavisPath)

    #Coverage
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_covVis.tsv")
    nPhaseFunctions.generateCoverageVis(simpleOutPath,datavisPath)

    #Discordance
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_discordanceVis.tsv")
    nPhaseFunctions.generateDiscordanceVis(simpleOutPath,datavisPath)

    readmeText="\nPlot can be found at "+datavisFolderPath
    print(readmeText)
    updateLog(readmePath,readmeText)


    #Generate FastQ Files
    haplotigReadNameFile=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_clusterReadNames.tsv")
    fastQFilePrefix=os.path.join(phasedFastqPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_")

    nPhaseFunctions.generateLongReadFastQFiles(haplotigReadNameFile,args.longReadFile,fastQFilePrefix)#Needs to work more efficiently, shouldn't just load the entire fastq into memory...right? (why not?) Do I do that at any other point? #Make this last just in case.

    readmeText="\nLong reads can be found in "+phasedFastqPath
    print(readmeText)
    updateLog(readmePath,readmeText)

    print("You can consult the readme at "+readmePath+" if you want a bit of guidance about your results. Please raise any issues on https://github.com/nPhasePipeline/nPhase")

    return 0

    ####################


def nPhaseAlgorithm(args):

    basePath=os.path.join(args.outputFolder,args.strainName)
    phasedPath=os.path.join(basePath,"Phased")
    phasedFastqPath=os.path.join(basePath,"Phased","FastQ")
    datavisFolderPath=os.path.join(basePath,"Phased","Plots")
    logPath=os.path.join(basePath,"Logs") #Currently not producing any logs

    allPaths=[basePath,phasedPath,phasedFastqPath,datavisFolderPath,logPath]

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    readmePath=os.path.join(basePath,"Readme.txt")
    fullLogPath=os.path.join(basePath,"Logs","fullLog.txt")

    #Run everything through nPhase
    phaseTool.nPhase(args.validatedSNPAssignmentsFile,args.strainName,args.contextDepthsFile,phasedPath,basePath,args.reference,args.minSim,args.minOvl,args.minLen,args.maxID)

    readmeText="\nPhased files can be found at "+phasedPath+"\nThe *_variants.tsv file contains information on the consensus heterozygous variants present in each predicted haplotig.\nThe *_clusterReadNames.tsv file contains information on the reads which comprise each cluster."
    print(readmeText)
    updateLog(readmePath,readmeText)

    #Simplify datavis
    dataVisPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_phasedDataFull.tsv")
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_phasedDataSimple.tsv")
    nPhaseFunctions.simplifyDataVis(dataVisPath,simpleOutPath,1000)

    #Generate plots
    datavisPath=os.path.join(datavisFolderPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_")
    nPhaseFunctions.generatePhasingVis(simpleOutPath,datavisPath)

    #Coverage
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_covVis.tsv")
    nPhaseFunctions.generateCoverageVis(simpleOutPath,datavisPath)

    #Discordance
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_discordanceVis.tsv")
    nPhaseFunctions.generateDiscordanceVis(simpleOutPath,datavisPath)

    readmeText="\nPlot can be found at "+datavisFolderPath
    print(readmeText)
    updateLog(readmePath,readmeText)

    #Generate FastQ Files
    haplotigReadNameFile=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_clusterReadNames.tsv")
    fastQFilePrefix=os.path.join(phasedFastqPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_")

    nPhaseFunctions.generateLongReadFastQFiles(haplotigReadNameFile,args.longReadFile,fastQFilePrefix) #Needs to work more efficiently, shouldn't just load the entire fastq into memory...right? (why not?) Do I do that at any other point? #Make this last just in case.

    readmeText="\nLong reads can be found in "+phasedFastqPath
    print(readmeText)
    updateLog(readmePath,readmeText)

    print("You can consult the readme at "+readmePath+" if you want a bit of guidance about your results. Please raise any issues on https://github.com/nPhasePipeline/nPhase")

    return 0

def partialPipeline(args):
    allPaths=[]

    #Create folder structure
    basePath=os.path.join(args.outputFolder,args.strainName)
    mappedShortReadPath=os.path.join(basePath,"Mapped","shortReads")
    if args.mappedLongReads=="noMapped":
        mappedLongReadPath=os.path.join(basePath,"Mapped","longReads")
        allPaths.append(mappedLongReadPath)
    variantCalledShortReadPath=os.path.join(basePath,"VariantCalls","shortReads")
    variantCalledLongReadPath=os.path.join(basePath,"VariantCalls","longReads")
    overlapPath=os.path.join(basePath,"Overlaps")
    phasedPath=os.path.join(basePath,"Phased")
    phasedFastqPath=os.path.join(basePath,"Phased","FastQ")
    datavisFolderPath=os.path.join(basePath,"Phased","Plots")
    logPath=os.path.join(basePath,"Logs") #Currently not producing any logs

    allPaths+=[basePath,mappedShortReadPath,variantCalledShortReadPath,variantCalledLongReadPath,overlapPath,phasedPath,phasedFastqPath,datavisFolderPath,logPath]

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    readmePath=os.path.join(basePath,"Readme.txt")
    fullLogPath=os.path.join(basePath,"Logs","fullLog.txt")

    ########################
    #Pre-process reference #
    ########################

    #Make sure the reference is indexed (might be missing a process)
    if not os.path.isfile(args.reference+".fai"):
        p=subprocess.run(["samtools","faidx",args.reference],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["samtools","faidx",args.reference])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)

    if not os.path.isfile(args.reference+".bwt"):
        p=subprocess.run(["bwa","index",args.reference],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["bwa","index",args.reference])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)

    if not os.path.isfile(args.reference+".dict"):
        p=subprocess.run(["gatk","CreateSequenceDictionary","-R",args.reference],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["gatk","CreateSequenceDictionary","-R",args.reference])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)

    ########################
    #Pre-process long reads#
    ########################

    #Map long reads to reference
    if args.mappedLongReads=="noMapped":
        splitReadFlag="260" #This flag allows split reads
        outputLog,systemMessage=nPhaseFunctions.longReadMapping(args.strainName,args.longReadFile,args.reference,mappedLongReadPath,splitReadFlag,args.longReadPlatform,args.threads)
        print(systemMessage)
        updateLog(fullLogPath,outputLog)


    #########################
    #Pre-process short reads#
    #########################

    #Map short reads to reference
    if args.mappedShortReads=="noMapped" and args.vcfFile=="noVCF":
        outputLog,systemMessage=nPhaseFunctions.shortReadMapping(args.strainName,args.shortReadFile_R1,args.shortReadFile_R2,args.reference,mappedShortReadPath,args.threads)
        print(systemMessage)
        updateLog(fullLogPath,outputLog)

    if args.vcfFile=="noVCF":
        #Variant call short reads and select SNPs only
        estimatedPloidy="2" #This argument is required for GATK's variant calling and we can expect most SNPs to only have two alleles anyway

        if args.mappedShortReads=="noMapped":
            shortReadBam=os.path.join(mappedShortReadPath,args.strainName+".final.bam")
        else:
            shortReadBam=args.mappedShortReads
        shortReadVCF=os.path.join(variantCalledShortReadPath,args.strainName+".vcf")
        p=subprocess.run(["gatk","HaplotypeCaller","-R",args.reference,"-ploidy",estimatedPloidy,"-I",shortReadBam,"-O",shortReadVCF],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
        logText="COMMAND: "+" ".join(["gatk","HaplotypeCaller","-R",args.reference,"-ploidy",estimatedPloidy,"-I",shortReadBam,"-O",shortReadVCF])+"\n\n"
        if p.stderr!="" or p.stdout !="":
            logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
        updateLog(fullLogPath,logText)
    else:
        shortReadVCF=args.vcfFile

    shortReadSNPsVCF=os.path.join(variantCalledShortReadPath,args.strainName+".SNPs.vcf")
    p=subprocess.run(["gatk","SelectVariants","-R",args.reference,"--variant",shortReadVCF,"-O",shortReadSNPsVCF,"--select-type-to-include","SNP"],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    logText="COMMAND: "+" ".join(["gatk","SelectVariants","-R",args.reference,"--variant",shortReadVCF,"-O",shortReadSNPsVCF,"--select-type-to-include","SNP"])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
    updateLog(fullLogPath,logText)


    #Extract heterozygous positions from VCF file
    #Ugh this is a whole thing.

    SNPVCFFilePath=os.path.join(variantCalledShortReadPath,args.strainName+".SNPs.vcf")
    hetSNPVCFFilePath=os.path.join(variantCalledShortReadPath,args.strainName+".hetSNPs.vcf")

    hetSNPVCFText=""

    SNPVCFFile=open(SNPVCFFilePath,"r")
    for line in SNPVCFFile:
        if "#" in line:
            hetSNPVCFText+=line
        elif "AF=1.00" not in line:
            line=line.strip("\n")
            line=line.replace(";","\t")
            line=line.split("\t")
            line=[line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7]]
            hetSNPVCFText+="\t".join(line)+"\n"

    SNPVCFFile.close()

    hetSNPVCFFile=open(hetSNPVCFFilePath,"w")
    hetSNPVCFFile.write(hetSNPVCFText)
    hetSNPVCFText="" #Just saving some memory since this isn't in a function
    hetSNPVCFFile.close()

    print("Identified heterozygous SNPs in short read VCF")

    shortReadPositionsOutputFilePath=os.path.join(mappedShortReadPath,args.strainName+".hetSNPs.positions.tsv")
    nPhaseFunctions.getShortReadPositions(hetSNPVCFFilePath,shortReadPositionsOutputFilePath)

    shortReadSNPsBedFilePath=os.path.join(mappedShortReadPath,args.strainName+".hetSNPs.bed")

    shortReadSNPsBedText=""

    shortReadPositionsOutputFile=open(shortReadPositionsOutputFilePath,"r")
    for line in shortReadPositionsOutputFile:
        if "#" in line:
            shortReadSNPsBedText+=line
        elif "AF=1.00" not in line:
            line=line.strip("\n")
            line=line.replace(":","\t")
            line=line.split("\t")
            line=[line[0],str(int(line[1])-1),line[1]]
            shortReadSNPsBedText+="\t".join(line)+"\n"

    shortReadPositionsOutputFile.close()

    shortReadSNPsBedFile=open(shortReadSNPsBedFilePath,"w")
    shortReadSNPsBedFile.write(shortReadSNPsBedText)
    shortReadSNPsBedText="" #Just saving some memory since this isn't in a function
    shortReadSNPsBedFile.close()

    print("Extracted heterozygous SNP info based on short read VCF")

    #Reduce long reads to their heterozygous SNPs
    if args.mappedLongReads=="noMapped":
        cleanLongReadSamFile=os.path.join(mappedLongReadPath,args.strainName+".sorted.sam")
    else:
        cleanLongReadSamFile=args.mappedLongReads
    minQ=0.01   # Currently
    minMQ=0     # Not
    minAln=0.5  # Used
    longReadPositionNTFile=os.path.join(variantCalledLongReadPath,args.strainName+".hetPositions.SNPxLongReads.tsv")
    nPhaseFunctions.assignLongReadToSNPs(cleanLongReadSamFile,shortReadSNPsBedFilePath,args.reference,minQ,minMQ,minAln,longReadPositionNTFile)

    #Only keep the longest reads and get rid of duplicate heterozygous SNP profiles, compensate by keeping context coverage information.
    minCov=0     # Currently
    minRatio=0   # Not
    minTrioCov=0 # Used
    validatedSNPAssignmentsFile=os.path.join(variantCalledLongReadPath,args.strainName+".hetPositions.SNPxLongReads.validated.tsv")
    contextDepthsFile=os.path.join(overlapPath,args.strainName+".contextDepths.tsv")
    nPhaseFunctions.longReadValidation(longReadPositionNTFile,minCov,minRatio,minTrioCov,validatedSNPAssignmentsFile,contextDepthsFile)

    #Run everything through nPhase
    phaseTool.nPhase(validatedSNPAssignmentsFile,args.strainName,contextDepthsFile,phasedPath,basePath,args.reference,args.minSim,args.minOvl,args.minLen,args.maxID)

    readmeText="\nPhased files can be found at "+phasedPath+"\nThe *_variants.tsv file contains information on the consensus heterozygous variants present in each predicted haplotig.\nThe *_clusterReadNames.tsv file contains information on the reads which comprise each cluster."
    print(readmeText)
    updateLog(readmePath,readmeText)

    #Simplify datavis
    dataVisPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_phasedDataFull.tsv")
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_phasedDataSimple.tsv")
    nPhaseFunctions.simplifyDataVis(dataVisPath,simpleOutPath,1000)

    ################
    #Generate plots#
    ################

    #Phased
    datavisPath=os.path.join(datavisFolderPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_")
    nPhaseFunctions.generatePhasingVis(simpleOutPath,datavisPath)

    #Coverage
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_covVis.tsv")
    nPhaseFunctions.generateCoverageVis(simpleOutPath,datavisPath)

    #Discordance
    simpleOutPath=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_discordanceVis.tsv")
    nPhaseFunctions.generateDiscordanceVis(simpleOutPath,datavisPath)

    readmeText="\nPlot can be found at "+datavisFolderPath
    print(readmeText)
    updateLog(readmePath,readmeText)


    #Generate FastQ Files
    haplotigReadNameFile=os.path.join(phasedPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_clusterReadNames.tsv")
    fastQFilePrefix=os.path.join(phasedFastqPath,args.strainName+"_"+str(args.minOvl)+"_"+str(args.minSim)+"_"+str(args.maxID)+"_"+str(args.minLen)+"_")

    nPhaseFunctions.generateLongReadFastQFiles(haplotigReadNameFile,args.longReadFile,fastQFilePrefix)#Needs to work more efficiently, shouldn't just load the entire fastq into memory...right? (why not?) Do I do that at any other point? #Make this last just in case.

    readmeText="\nLong reads can be found in "+phasedFastqPath
    print(readmeText)
    updateLog(readmePath,readmeText)

    print("You can consult the readme at "+readmePath+" if you want a bit of guidance about your results. Please raise any issues on https://github.com/nPhasePipeline/nPhase")

    return 0

def cleaning(args):
    cleaningTool.nPhaseCleaning(args.phasingResult,args.longReadFile,args.strainName,args.percentKept,args.maxDiscordance,args.deduplicate,args.FFBool)

def main():
    parser=argparse.ArgumentParser(prog="nPhase",description='Full ploidy agnostic phasing pipeline',add_help=False)
    parser.add_argument('--version',action='version',version='%(prog)s 1.1.3')

    #CREATING SEPARATE MODES
    subparsers = parser.add_subparsers(help="selecting 'pipeline' will run through all of the steps, 'nPhase' will only perform the phasing operation, 'partial' will only perform part of the pipeline, 'cleaning' will run automated cleaning steps on nPhase results",dest='command')

    parser_a = subparsers.add_parser('pipeline', help='Run the entire nPhase pipeline on your sample.',parents=[parser])
    required_a=parser_a.add_argument_group('required arguments')

    parser_b = subparsers.add_parser('algorithm', help='Only run the nPhase algorithm. NOTE: This will require files generated by running the pipeline mode',parents=[parser])
    required_b=parser_b.add_argument_group('required arguments')

    parser_c = subparsers.add_parser('partial', help='Only run parts of the nPhase pipeline on your sample. NOTE: if you run into any issues with this mode please use the pipeline mode instead.',parents=[parser])
    required_c=parser_c.add_argument_group('required arguments')

    cleaningParser=subparsers.add_parser('cleaning',help='Clean nPhase results using a simple three step pipeline. Filtering, merging, gap filling.',parents=[parser])
    required_Cleaning=cleaningParser.add_argument_group('required arguments')

    #SHARED REQUIRED ARGUMENTS

    for parserReference in [required_a,required_b,required_c,required_Cleaning]:
        parserReference.add_argument('--sampleName',required=True, dest='strainName',help='Name of your sample, ex: "Individual_1"')
        parserReference.add_argument('--longReads',required=True, dest='longReadFile',help='Path to long read FastQ file, ex: /home/longReads/Individual_1.fastq.gz')

    #FULL PIPELINE SPECIFIC ARGUMENTS
    required_a.add_argument('--longReadPlatform',required=True, dest='longReadPlatform',choices=['ont','pacbio'],help="Long read platform, must be 'ont' or 'pacbio'")
    required_a.add_argument('--R1',required=True, dest='shortReadFile_R1',help='Path to paired end short read FastQ file #1, ex: /home/shortReads/Individual_1_R1.fastq.gz')
    required_a.add_argument('--R2',required=True, dest='shortReadFile_R2',help='Path to paired end short read FastQ file #2, ex: /home/shortReads/Individual_1_R2.fastq.gz')

    #NPHASE SPECIFIC ARGUMENTS
    required_b.add_argument('--contextDepth',required=True, dest='contextDepthsFile',help='Path to context depths file, ex: /home/phased/Individual_1/Overlaps/Individual_1.contextDepths.tsv')
    required_b.add_argument('--processedLongReads',required=True, dest='validatedSNPAssignmentsFile',help='Path to validated long read SNPs, ex: /home/phased/Individual_1/VariantCalls/longReads/Individual_1.hetPositions.SNPxLongReads.validated.tsv')

    #PARTIAL PIPELINE SPECIFIC ARGUMENTS
    required_c.add_argument('--mappedShortReads', dest='mappedShortReads', default="noMapped",help='Path to mapped, sorted short read file in BAM format, ex: /home/phased/Individual_1/Mapped/shortReads/Individual_1.final.bam')
    required_c.add_argument('--vcf',dest='vcfFile',default="noVCF",help="Path to VCF file (must be generated with GATK --ploidy=2), ex:/home/phased/Individual_1/VariantCalls/shortReads/Individual.vcf")
    required_c.add_argument('--mappedLongReads', dest='mappedLongReads', default="noMapped",help='Path to mapped long read file in SAM format and should be filtered to keep split reads (flag 260), ex: /home/phased/Individual_1/Mapped/longReads/Individual_1.sorted.sam')
    required_c.add_argument('--longReadPlatform', dest='longReadPlatform',default="noLRMapping",choices=['ont','pacbio'],help="Long read platform, must be 'ont' or 'pacbio'")
    required_c.add_argument('--R1', dest='shortReadFile_R1',default="noSRMapping",help='Path to paired end short read FastQ file #1, ex: /home/shortReads/Individual_1_R1.fastq.gz')
    required_c.add_argument('--R2', dest='shortReadFile_R2',default="noSRMapping",help='Path to paired end short read FastQ file #2, ex: /home/shortReads/Individual_1_R2.fastq.gz')

    #CLEANING SPECIFIC ARGUMENTS
    required_Cleaning.add_argument('--resultFolder',required=True,dest='phasingResult',help='Path to output folder, ex: /home/nphaseResults/Individual_1/')
    required_Cleaning.add_argument('--filterPct', dest='percentKept',type=float,nargs="?",default=2.0,help='Percentage of results to filter out. Default 2.0')
    required_Cleaning.add_argument('--fillGaps',dest='deduplicate',type=bool,nargs="?",default=True,help='Attempt to fill gaps by redistributing reads from most covered cluster. Set to 0 to disable this step. Default 1')
    required_Cleaning.add_argument('--filterFirst',dest='FFBool',type=bool,nargs="?",default=False,help='Perform the filtering step before the merging step. Set to 1 to enable. Default 0')
    required_Cleaning.add_argument('--setMaxDiscordance',dest='maxDiscordance',type=float,nargs="?",default=0.0,help='If this is set higher than 0, use this percentage as the stopping rule for merging instead of using the mean discordance of raw results. Inputting 6.8 means 6.8%% discordance max allowed. Default 0')

    #OPTIONAL SHARED ARGUMENTS
    for parserReference in [required_a,required_b,required_c]:
        parserReference.add_argument('--reference',required=True,dest='reference',help='Path to fasta file of reference genome to align to, ex: /home/reference/Individual_reference.fasta')
        parserReference.add_argument('--output',required=True,dest='outputFolder',help='Path to output folder, ex: /home/phased/')
        parserReference.add_argument('--threads',dest='threads',type=str,nargs="?",default="8",help='Number of threads to use on some steps. Default 8')
        parserReference.add_argument('--maxID',dest='maxID',type=float,nargs="?",default=0.05,help='MaxID parameter, determines how different two clusters must be to prevent them from merging. Default 0.05')
        parserReference.add_argument('--minOvl',dest='minOvl',type=float,nargs="?",default=0.1,help='minOvl parameter, determines the minimal percentage of overlap required to allow a merge between two clusters that have fewer than 100 heterozygous SNPs in common. Default 0.1')
        parserReference.add_argument('--minSim',dest='minSim',type=float,nargs="?",default=0.01,help='minSim parameter, determines the minimal percentage of similarity required to allow a merge between two clusters. Default 0.01')
        parserReference.add_argument('--minLen',dest='minLen',type=int,nargs="?",default=0,help='minLen parameter, any cluster based on fewer than N reads will not be output. Default 0')

    args, unknown = parser.parse_known_args()

    if "mappedShortReads" in dir(args) or "vcf" in dir(args) or "mappedLongReads" in dir(args):
        partialPipeline(args)
    elif "longReadPlatform" in dir(args):
        nPhasePipeline(args)
    elif "contextDepthsFile" in dir(args):
        nPhaseAlgorithm(args)
    elif "phasingResult" in dir(args):
        cleaning(args)
    else:
        print("Select which mode of nphase to use. For example:\n\nnphase pipeline [args] [-h]\n\nnphase algorithm [args] [-h]\n\nnphase partial [args] [-h]\n\nnphase filter [args] [-h]")

    exit(0)

if __name__=="__main__":
    main()
    exit(0)



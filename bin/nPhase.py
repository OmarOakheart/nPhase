import operator
import os
from itertools import combinations
import sys
import sortedcontainers

def identity(readIDict,clusterAID,clusterBID,commonPositions):
    allBaseDict={}
    for position in commonPositions:
        baseA=clusterAID[position]
        baseB=clusterBID[position]
        basePossibilities=set().union(baseA["demo"].keys(),baseB["demo"].keys())
        baseNew={"demo":{},"stats":{"N":baseA["stats"]["N"]+baseB["stats"]["N"],"total":baseA["stats"]["total"]+baseB["stats"]["total"]}}
        for base in basePossibilities:
            if baseNew["stats"]["total"]>0:
                baseNew["demo"][base]=(baseA["demo"].get(base,0)*baseA["stats"]["total"]+baseB["demo"].get(base,0)*baseB["stats"]["total"])/baseNew["stats"]["total"]
            else:
                baseNew["demo"][base]=0
        allBaseDict[position]=baseNew
    for position in clusterAID.keys():
        if position not in commonPositions:
            allBaseDict[position]=clusterAID[position]
    for position in clusterBID.keys():
        if position not in commonPositions:
            allBaseDict[position]=clusterBID[position]
    return allBaseDict

def identityChangeBool(clusterID,mergedClusterID,commonPositions,maxID): #Penalize both
    overallChange=0
    allSeqNumbers=0
    if len(commonPositions)==0: #I don't like having this here, but necessary to hack the final stitching, to clean up later. Is it still necessary?
        return False
    for position in commonPositions:
        demographics=clusterID[position]
        bestDemo=max(demographics["demo"].values())
        allSeqNumbers+=mergedClusterID[position]["stats"]["N"]
        for demo, proportion in demographics["demo"].items():
            if proportion==bestDemo:
                if mergedClusterID[position]["demo"][demo]-proportion<0: #If it's the best, and after the merge the percentage decreased, then we want to consider that a change in terms of the identity of the cluster
                    overallChange+=abs(mergedClusterID[position]["demo"][demo]-proportion)*mergedClusterID[position]["stats"]["N"]
            else:
                if mergedClusterID[position]["demo"][demo]-proportion>0: #If it isn't the best, and after the merge the percentage increased, then we want to consider that a change in terms of the identity of the cluster
                    overallChange+=abs(proportion-mergedClusterID[position]["demo"][demo])*mergedClusterID[position]["stats"]["N"]
    if allSeqNumbers==0: #Disallow the rare case that two clusters have something in common that isn't supported by any contextDepth info
        return True
    if overallChange/allSeqNumbers>maxID:
        return True
    else:
        return False

def identityTheftBool(positionsA,positionsB,clusterA,clusterB,consensusA,consensusB,clusterAID,clusterBID,contextDepths,readIDict,maxID): #There are two possible ways to do this, either by discouraging dissidence, or deposition
    #This is the dissidence mode. Deposition mode would take into account only successful actions against established bases (much more aggressive clustering).
    commonPositions=positionsA&positionsB
    mergedClusterID=identity(readIDict,clusterAID,clusterBID,commonPositions)
    AChanges=identityChangeBool(clusterAID,mergedClusterID,commonPositions,maxID)
    if identityChangeBool(clusterAID,mergedClusterID,commonPositions,maxID) and identityChangeBool(clusterBID,mergedClusterID,commonPositions,maxID):
        return True
    else:
        return False

def consensus(clusterID,contextDepths):
    consensus=set()
    for position in clusterID.keys():
        highestPct=max(clusterID[position]["demo"].values())
        for base, percentage in clusterID[position]["demo"].items():
            if percentage==highestPct: #To ensure that if there's a 50/50 split, they both get a chance.
                consensus.add(position+"="+base)
    return frozenset(consensus)

def consensusInitializer(sequenceList,contextDepths):
    basePositions=set()
    for sequence in sequenceList:
        basePositions.update([x.split("=")[0] for x in sequence])
    baseDepths={}
    for basePos in basePositions:
        baseDepths[basePos]={}
    for sequence in sequenceList: #Each sequence votes once, therefore the same context is present multiple times in the same cluster results in a stronger vote, hopefully enough to overwhelm chimera but too little to amplify errors
        for base in sequence:
            basePos,baseNt=base.split("=")
            for context, contextCount in contextDepths[base].items():
                if context.issubset(sequence):
                    baseDepths[basePos][baseNt]=baseDepths[basePos].setdefault(baseNt,0)+contextCount
    consensus=set()
    for base, baseDepth in baseDepths.items():
        best=max(baseDepth.values(),default=0)
        for baseNt, depth in baseDepth.items():
            if depth==best:
                consensus.add(base+"="+baseNt)
    return frozenset(consensus)

def bestTwo(sequences,similarityIndex,minSim,contextDepths,maxID,readIDict,bannedClusterNames):
    bestCouple=[]
    badCombinations=[]
    maxLen=0
    indx=0
    indexes=[]
    for combination in reversed(similarityIndex):
        indx-=1
        nameA=combination[2]
        nameB=combination[3]
        if nameA in bannedClusterNames or nameB in bannedClusterNames:
            indexes.append(indx)
        else:
            cacheA=sequences[combination[2]]
            cacheB=sequences[combination[3]]
            similarity=combination[0]
            if similarity>=minSim:
                if not identityTheftBool(cacheA["positions"],cacheB["positions"],cacheA["names"],cacheB["names"],cacheA["consensus"],cacheB["consensus"],cacheA["clusterID"],cacheB["clusterID"],contextDepths,readIDict,maxID):
                    return (combination[2],combination[3]),indexes
                else:
                    indexes.append(indx)
            elif similarity<minSim:
                break
    return [],indexes

def getSimilarity(consensusA,consensusB,commonPos,minOvl,minSim):
    if len(consensusA)<len(consensusB):
        setA=consensusA
        setB=consensusB
    else:
        setA=consensusB
        setB=consensusA
    if len(commonPos)/len(setA)>=minOvl or len(commonPos)>100:
        commonSet=setA&setB
        localSimilarity=len(commonSet)/len(commonPos)
        if localSimilarity<minSim:
            localSimilarity=0
    else:
        localSimilarity=0
    return localSimilarity

def initializeCache(readIDict,name,contextDepths):
    cache={}
    i=0
    for cluster in readIDict.keys():
        cluster=[cluster]
        seqNames=set()
        seqCluster=set()
        for readName in cluster:
            sequence=frozenset([k+"="+next(iter(v["demo"])) for k,v in readIDict[readName].items()])
            seqNames.add(readName)
            seqCluster.add(sequence)
        consensusSeq=consensusInitializer(seqCluster,contextDepths)
        SNPPositions=set([SNP.split("=")[0] for SNP in consensusSeq])
        if len(seqNames)==1:
            clusterID=readIDict[list(seqNames)[0]]
        else:
            print("This broke")
            exit(1)
        cache[name+str(i)]={"names":seqNames,"cluster":seqCluster,"clusterID":clusterID,"consensus":consensusSeq,"positions":SNPPositions,"similarities":{},"overlaps":set()}
        i+=1
    return cache

def fillCache(cache,minOvl,minSim):
    for combination in combinations(list(cache.keys()),2):
        cacheA=cache[combination[0]]
        cacheB=cache[combination[1]]
        commonPos=cacheA["positions"]&cacheB["positions"]
        if len(commonPos)>0:
            cacheA["overlaps"].add(combination[1])
            cacheB["overlaps"].add(combination[0])
            similarity=getSimilarity(cacheA["consensus"],cacheB["consensus"],commonPos,minOvl,minSim)
            if similarity>=minSim:
                if len(cacheA["consensus"])<len(cacheB["consensus"]):
                    cacheA["similarities"][combination[1]]=similarity
                else:
                    cacheB["similarities"][combination[0]]=similarity
    return cache

def generateSimilarityIndex(cache):
    similarityIndex=sortedcontainers.SortedList()
    for simCluster in cache.keys():
        for combination, score in cache[simCluster]["similarities"].items():
            similarPositions=cache[simCluster]["positions"]&cache[combination]["positions"]
            similarityIndex.add([score,len(similarPositions),simCluster,combination])
    return similarityIndex

def clusterInterlockChimera(minOvl,readIDict,maxID,minLen,minSim,contextDepths):
    prevLen=0
    allReadIDict={}
    allReadIDict.update(readIDict)
    #Initializing cachedCluster
    print("Initializing cachedCluster")
    cachedSimilarities=initializeCache(readIDict,"cluster",contextDepths)
    #Filling cachedCluster with similarity info
    print("Filling cachedCluster with similarity information")
    cachedSimilarities=fillCache(cachedSimilarities,minOvl,minSim)
    print("Preparing initial similarity index")
    bannedClusterNames=set()
    similarityIndex=generateSimilarityIndex(cachedSimilarities)
    print("Starting clustering loop ("+str(len(cachedSimilarities.keys()))+" sequences)")
    chrSizes=len(cachedSimilarities.keys())
    i=0
    while len(cachedSimilarities.keys())!=prevLen:
        prevLen=len(cachedSimilarities.keys())
        pair,indexes=bestTwo(cachedSimilarities,similarityIndex,minSim,contextDepths,maxID,allReadIDict,bannedClusterNames)
        bannedClusterNames.update(pair)
        for indx in reversed(indexes):
            similarityIndex.pop(indx)
        if pair!=[]:
            first=pair[0]
            second=pair[1]
            newName="mergedCluster_"+str(i)
            newNames=cachedSimilarities[first]["names"]|cachedSimilarities[second]["names"]
            newCluster=cachedSimilarities[first]["cluster"]|cachedSimilarities[second]["cluster"]
            commonPositions=set([x.split("=")[0] for x in cachedSimilarities[first]["consensus"]])&set([x.split("=")[0] for x in cachedSimilarities[second]["consensus"]])
            newClusterID=identity(allReadIDict,cachedSimilarities[first]["clusterID"],cachedSimilarities[second]["clusterID"],commonPositions)
            newOverlaps=cachedSimilarities[first]["overlaps"]|cachedSimilarities[second]["overlaps"]
            newOverlaps.remove(first)
            newOverlaps.remove(second)
            for overlap in newOverlaps:
                cachedSimilarities[overlap]["overlaps"].add(newName)
            newConsensus=consensus(newClusterID,contextDepths)
            newSNPPositions=set([SNP.split("=")[0] for SNP in newConsensus])
            newCache={"names":newNames,"cluster":newCluster,"clusterID":newClusterID,"consensus":newConsensus,"positions":newSNPPositions,"similarities":{},"overlaps":newOverlaps}
            del cachedSimilarities[first]
            del cachedSimilarities[second]
            for cache in newCache["overlaps"]:
                cacheA=cachedSimilarities[cache]
                cacheB=newCache
                if first in cacheA["similarities"]:
                    del cacheA["similarities"][first]
                if second in cacheA["similarities"]:
                    del cacheA["similarities"][second]
                if first in cacheA["overlaps"]:
                    cacheA["overlaps"].remove(first)
                if second in cacheA["overlaps"]:
                    cacheA["overlaps"].remove(second)
                newCommonPos=cacheA["positions"]&cacheB["positions"]
                if len(newCommonPos)>minOvl:
                    similarity=getSimilarity(cacheA["consensus"],cacheB["consensus"],newCommonPos,minOvl,minSim)
                else:
                    similarity=0
                if similarity>=minSim:
                    if len(cacheA["consensus"])<len(cacheB["consensus"]):
                        cacheA["similarities"][newName]=similarity
                        similarityIndex.add([similarity,len(newCommonPos),cache,newName])
                    else:
                        cacheB["similarities"][cache]=similarity
                        similarityIndex.add([similarity,len(newCommonPos),newName,cache])
            cachedSimilarities[newName]=newCache
            i+=1
    #This is where you get the results output, needs another shell script or a python script to do the rest
    actualClusters={}
    resultClusters={}
    for cacheName, cacheInfo in cachedSimilarities.items():
        if len(cacheInfo["names"])>minLen: #Only allow clusters with at least this many reads
            actualClusters[cacheName]=cacheInfo["consensus"]
            resultClusters[cacheName]=cacheInfo
    if len(actualClusters)>0:
        print(len(actualClusters),"clusters")
    else:
        print("No clusters.")
    return actualClusters,resultClusters

#############################
#SPLIT READ CODE STARTS HERE#
#############################

#Using the base cluster info to determine the N SNPs at the edges of the formed clusters

def getSNPPositions(cluster):
    SNPPositions=[]
    for SNP in cluster:
        SNPPositions.append(int(SNP.split("=")[0].split(":")[1]))
    SNPPositions.sort()
    contig=SNP.split(":")[0]
    sortedSNPPositions=[contig+":"+str(SNP) for SNP in SNPPositions]
    return sortedSNPPositions

def getBaseClusterEdges(N,clusters):
    baseClusterEdges=set()
    for cluster in clusters:
        sortedSNPPositions=getSNPPositions(cluster)
        startEdge=frozenset(sortedSNPPositions[0:min(N,len(cluster))])
        endEdge=frozenset(sortedSNPPositions[max(len(cluster)-N,0):len(cluster)])
        baseClusterEdges.add(startEdge)
        baseClusterEdges.add(endEdge)
    return baseClusterEdges

def mergeClusterEdges(clusterEdges):
    newClusterEdges=set()
    newClusterEdges.update(clusterEdges)
    used=set()
    change=True
    while change==True:
        change=False
        prevLen=len(newClusterEdges)
        for combination in combinations(newClusterEdges,2):
            clusterEdgeA=combination[0]
            clusterEdgeB=combination[1]
            if clusterEdgeA&clusterEdgeB!=set():
                change=True
                if clusterEdgeA in newClusterEdges:
                    newClusterEdges.remove(clusterEdgeA)
                if clusterEdgeB in newClusterEdges:
                    newClusterEdges.remove(clusterEdgeB)
                newClusterEdges.add(frozenset(clusterEdgeA|clusterEdgeB))
                break
    return newClusterEdges

def getSplitReadPositions(chimericReadIDict):
    splitReadPositions=set()
    for splitRead,positions in chimericReadIDict.items():
        splitReadPositions.add(frozenset(positions))
    return splitReadPositions

def getSplitEdgeOverlapProfiles(splitReads,mergedClusterEdges):
    splitReadProfiles=set()
    for splitReadPositions in splitReads:
        splitReadProfile=set()
        for edge in mergedClusterEdges:
            if len(edge&splitReadPositions)>0:
                splitReadProfile.add(edge)
        splitReadProfiles.add(frozenset(splitReadProfile))
    return splitReadProfiles

def getSplitReadProfile(splitReadPositions,mergedClusterEdges):
    splitReadProfile=set()
    for edge in mergedClusterEdges:
        if len(edge&splitReadPositions)>0:
            splitReadProfile.add(edge)
    return splitReadProfile

def flattenProfile(splitReadProfile):
    flattenedProfile=set()
    for edge in splitReadProfile:
        flattenedProfile.update(edge)
    return flattenedProfile

def keepUsefulSplitReadsChr(usefulSplitReadProfiles,chimericReadIDict,mergedClusterEdges):
    usefulSplitReads={}
    for splitRead,splitReadInfo in chimericReadIDict.items():
        splitReadPositions=set(splitReadInfo.keys())
        splitReadProfile=getSplitReadProfile(set(splitReadInfo.keys()),mergedClusterEdges)
        if splitReadProfile in usefulSplitReadProfiles:
            splitReadPositions=flattenProfile(splitReadProfile)
            newSplitReadInfo={}
            for position,info in splitReadInfo.items():
                if position in splitReadPositions:
                    newSplitReadInfo[position]=info
            chrms=set()
            for position in newSplitReadInfo:
                chrms.add(position.split(":")[0])
            chrmSNPs={}
            for chrm in chrms:
                chrmSNPs[chrm]={}
            for position, info in newSplitReadInfo.items():
                chrmSNPs[position.split(":")[0]][position]=info
            for chrm, posInfo in chrmSNPs.items():
                usefulSplitReads[splitRead]=posInfo
    return usefulSplitReads

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

def nPhase(longReadSNPAssignments,strainName,contextDepthsFilePath,outFolder,mainFolder,referenceFilePath,minSim,minOvl,minLen,maxID):

    #Loading files

    print("Loading reads.")

    SNPAssignmentFile=open(longReadSNPAssignments,"r")

    SNPAssignments={}

    for line in SNPAssignmentFile:
        line=line.strip("\n").split("\t")
        SNPAssignments[line[0]]=frozenset(line[1:])

    SNPAssignmentFile.close()

    print("Reads loaded.")

    #Checking if there are any split reads:
    for v in SNPAssignments.values():
        if len(set(SNP.split(":")[0] for SNP in v))>1:
            print("THERE ARE SPLIT READS IN:",v)
            break

    #######################
    #We need the context depth info:

    print("Loading context depth.")

    contextDepths={}

    contextDepthsLines=[]

    contextDepthsFile=open(contextDepthsFilePath,"r")
    for line in contextDepthsFile:
        line=line.strip("\n").split("\t")
        contextDepthsLines.append(line)

    contextDepthsFile.close()

    for line in contextDepthsLines:
        contextDepths[line[0]]={}
        for context, count in zip(line[1::2],line[2::2]):
            contextDepths[line[0]][frozenset(context.split(" "))]=int(count)

    print("Context depth loaded")

    #This dict will keep track of split reads that we're about to merge
    splitTrackingDict={}

    #Getting chimeric read names
    chimericNames=set()

    for readName, sequence in SNPAssignments.items():
        if len(readName.split("_VCSTTP_"))>1:
            baseName=readName.split("_VCSTTP_")[0]
            chimericNames.add(baseName)
            splitTrackingDict[baseName]={}

    print("Split reads identified.")

    #Applying the context depth info to all of the reads
    chimericReadIDict={}
    for readName, sequence in SNPAssignments.items():
        if len(readName.split("_VCSTTP_"))>1 or readName in chimericNames:
            if len(readName.split("_VCSTTP_"))>1:
                baseName=readName.split("_VCSTTP_")[0]
            if readName in chimericNames:
                baseName=readName
            allBaseDict={}
            if len(sequence)>=10:
                sequence=set(sequence)
                positions=[]
                for base in sequence:
                    chr=base.split(":")[0]
                    if chr not in splitTrackingDict[baseName].keys():
                        splitTrackingDict[baseName][chr]=set()
                    splitTrackingDict[baseName][chr].add(readName)
                    base=base.split("=")
                    positions.append(base[0])
                    allBaseDict.setdefault(base[0],{"demo":{base[1]:0},"stats":{"total":0,"N":0}}).update({"demo":{base[1]:0}})
                    for context, contextCount in contextDepths[base[0]+"="+base[1]].items():
                        if context.issubset(sequence):
                            allBaseDict[base[0]]["demo"][base[1]]+=contextCount
                            allBaseDict[base[0]]["stats"]["total"]+=contextCount
                            allBaseDict[base[0]]["stats"]["N"]+=1
                for basePos, bases in allBaseDict.items():
                    for base in allBaseDict[basePos]["demo"].keys():
                        allBaseDict[basePos]["demo"][base]=allBaseDict[basePos]["demo"][base]/allBaseDict[basePos]["stats"]["total"]
                chimericReadIDict.setdefault("chimeric-"+baseName,{}).update(allBaseDict)

    print("Split reads processed.")

    #Getting all reads
    readIDict={}
    for readName, sequence in SNPAssignments.items():
        if len(sequence)>=10:
            allBaseDict={}
            for base in sequence:
                base=base.split("=")
                allBaseDict.setdefault(base[0],{"demo":{base[1]:0},"stats":{"total":0,"N":0}}).update({"demo":{base[1]:0}})
                for context, contextCount in contextDepths[base[0]+"="+base[1]].items():
                    if context.issubset(sequence):
                        allBaseDict[base[0]]["demo"][base[1]]+=contextCount
                        allBaseDict[base[0]]["stats"]["total"]+=contextCount
                        allBaseDict[base[0]]["stats"]["N"]+=1
            for basePos, bases in allBaseDict.items():
                for base in allBaseDict[basePos]["demo"].keys():
                    allBaseDict[basePos]["demo"][base]=allBaseDict[basePos]["demo"][base]/allBaseDict[basePos]["stats"]["total"]
            readIDict[readName]=allBaseDict

    print("All reads processed.")

    baseClusters,fullClusters=clusterInterlockChimera(minOvl,readIDict,maxID,minLen,minSim,contextDepths)
    baseClusters=baseClusters.values()

    baseClusterEdges=getBaseClusterEdges(100,baseClusters)

    mergedClusterEdges=mergeClusterEdges(baseClusterEdges)

    splitReadPositions=getSplitReadPositions(chimericReadIDict)
    splitReadProfiles=getSplitEdgeOverlapProfiles(splitReadPositions,mergedClusterEdges)

    usefulSplitReadProfiles=set()
    for profile in splitReadProfiles:
        if len(profile)>=2:
            usefulSplitReadProfiles.add(profile)

    mergedSplitReadProfiles=mergeClusterEdges(usefulSplitReadProfiles)

    newChimericDict=keepUsefulSplitReadsChr(usefulSplitReadProfiles,chimericReadIDict,mergedClusterEdges)

    #print(newChimericDict["chimeric-CCN_ACA_BMB_0"])

    cleanAllReadIDict={}
    cleanAllReadIDict.update(readIDict)
    cleanAllReadIDict.update(newChimericDict)

    cleanAllClusters,cleanFullClusters=clusterInterlockChimera(minOvl,cleanAllReadIDict,maxID,minLen,minSim,contextDepths)

    #############
    #OUTPUT CODE#
    #############

    visDataTextFull=giveMeFullData(cleanAllClusters)
    visDataFilePath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_phasedDataFull.tsv")
    visDataFile=open(visDataFilePath,"w")
    visDataFile.write(visDataTextFull)
    visDataFile.close()

    #Create sets of fastQ read names, to be turned into FASTQ files by another script

    clusterReadText=""

    for clusterName, clusterInfo in cleanFullClusters.items():
        chrms=set()
        for SNP in clusterInfo["consensus"]:
            chrms.add(SNP.split(":")[0])
        if len(chrms)>2:
            print("MULTI-CHROMOSOME CLUSTER ERROR",chrms)
            exit(1)
        for readName in clusterInfo["names"]:
            readName=readName.replace("chimeric-","")
            if readName in splitTrackingDict:
                splitReads=splitTrackingDict[readName][list(chrms)[0]]
                for splitRead in splitReads:
                    clusterReadText+="\t".join([clusterName,splitRead])+"\n"
            else:
                clusterReadText+="\t".join([clusterName,readName])+"\n"

    clusterReadFilePath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_clusterReadNames.tsv")
    clusterReadFile=open(clusterReadFilePath,"w")
    clusterReadFile.write(clusterReadText)
    clusterReadFile.close()

    #A tab separated file that shows you which SNPs were predicted to be in the same haplotigs

    phasedSNPText=""

    for clusterName, clusterInfo in cleanFullClusters.items():
        for SNP in clusterInfo["consensus"]:
            chrom=SNP.split(":")[0]
            pos=SNP.split(":")[1].split("=")[0]
            base=SNP.split("=")[1]
            phasedSNPText+="\t".join([clusterName,chrom,pos,base])+"\n"

    phasedSNPFilePath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_variants.tsv")
    phasedSNPFile=open(phasedSNPFilePath,"w")
    phasedSNPFile.write(phasedSNPText)
    phasedSNPFile.close()

    #A tab separated file showing the coverage for each haplotig

    readFilePath=os.path.join(mainFolder,"VariantCalls","longReads",strainName+".hetPositions.SNPxLongReads.validated.tsv")
    outPath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_covVis.tsv")
    windowSize=10000
    minCov=generateCoverage(clusterReadFilePath,readFilePath,outPath,windowSize)

    outPath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_discordanceVis.tsv")
    windowSize=10000
    generateDiscordance(clusterReadFilePath,readFilePath,outPath)

    outPath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_minCov="+str(minCov)+"_filterVis.tsv")
    filterCoverage(clusterReadFilePath,readFilePath,outPath,minCov)

    return "Phasing over"



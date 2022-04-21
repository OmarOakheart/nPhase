import operator
import os
from itertools import combinations
import sys
import sortedcontainers
import bin.nPhasePipelineFunctions as nPhaseFunctions
import multiprocessing


def identity(clusterAID,clusterBID,commonPositions):
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
    for position in commonPositions: #Can be parallelized
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

def identityTheftBool(positionsA,positionsB,clusterAID,clusterBID,maxID): #There are two possible ways to do this, either by discouraging dissidence, or deposition
    #This is the dissidence mode. Deposition mode would take into account only successful actions against established bases (much more aggressive clustering).
    commonPositions=positionsA&positionsB
    mergedClusterID=identity(clusterAID,clusterBID,commonPositions)
    if identityChangeBool(clusterAID,mergedClusterID,commonPositions,maxID) and identityChangeBool(clusterBID,mergedClusterID,commonPositions,maxID):
        return True
    else:
        return False

def consensus(clusterID):
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

def nBestPairs(sequences,similarityIndex,minSim,maxID,bannedClusterNames,nPairs):
    indx=0
    indexes=[]
    nBest = []
    bannedOverlaps = set()
    for combination in reversed(similarityIndex):
        indx-=1
        nameA=combination[2]
        nameB=combination[3]
        if nameA in bannedClusterNames or nameB in bannedClusterNames:
            indexes.append(indx)
        else:
            if combination[2] not in bannedOverlaps and combination[3] not in bannedOverlaps:
                cacheA=sequences[combination[2]]
                cacheB=sequences[combination[3]]
                similarity=combination[0]
                if similarity>=minSim:
                    if not identityTheftBool(cacheA["positions"],cacheB["positions"],cacheA["clusterID"],cacheB["clusterID"],maxID):
                        nBest.append([combination[2], combination[3]])
                        bannedOverlaps.update(cacheA["overlaps"], cacheB["overlaps"])
                        if len(nBest)==nPairs:
                            return nBest, indexes
                    else:
                        indexes.append(indx)
                elif similarity<minSim:
                    pass
    return nBest,indexes

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

def generateSimilarityIndex(cache):
    similarityIndex=sortedcontainers.SortedList()
    for simCluster in cache.keys():
        for combination, score in cache[simCluster]["similarities"].items():
            similarPositions=cache[simCluster]["positions"]&cache[combination]["positions"]
            similarityIndex.add([score,len(similarPositions),simCluster,combination])
    return similarityIndex

def cacheFiller(queue,outQueue,minOvl,minSim):
    localCache={}
    while True:
        combinationBatch=queue.get()
        if combinationBatch == "end":
            outQueue.put(localCache)
            return
        for combination in combinationBatch:
            localCache.setdefault(combination[0], {"overlaps": set(), "similarities": {}})
            localCache.setdefault(combination[1], {"overlaps": set(), "similarities": {}})
            cacheKey = combination[0]
            otherKey = combination[1]
            cacheA = cachedSimilarities[cacheKey]
            cacheB = cachedSimilarities[otherKey]
            commonPos = cacheA["positions"] & cacheB["positions"]
            if len(commonPos) > 0:
                localCache[cacheKey]["overlaps"].add(otherKey)
                localCache[otherKey]["overlaps"].add(cacheKey)
                similarity = getSimilarity(cacheA["consensus"], cacheB["consensus"], commonPos, minOvl, minSim)
                if similarity >= minSim:
                    if len(cacheA["consensus"]) < len(cacheB["consensus"]):
                        localCache[cacheKey]["similarities"][otherKey] = similarity
                    else:
                        localCache[otherKey]["similarities"][cacheKey] = similarity

def startCacheFillers(queue, outQueue, cachedSimilarityQueue, cachedSimilarities, minOvl, minSim, threadNumber):
    allProcesses = list()
    for thread in range(0, threadNumber-1):
        cacheFiller_p = multiprocessing.Process(target=cacheFiller, args=((queue),(outQueue),(minOvl),(minSim),))
        cacheFiller_p.daemon = True
        cacheFiller_p.start()
        allProcesses.append(cacheFiller_p)
    cacheMerger_p=multiprocessing.Process(target=cacheMerger, args=((outQueue),(cachedSimilarityQueue),(cachedSimilarities),))
    cacheMerger_p.daemon = True
    cacheMerger_p.start()
    allProcesses.append(cacheMerger_p)

    return allProcesses

def cacheMerger(outQueue,cachedSimilarityQueue,cachedSimilarities):
    while True:
        localDict = outQueue.get()
        if localDict=="end":
            cachedSimilarityQueue.put(cachedSimilarities)
            return
        for cacheName in localDict.keys():
            cachedSimilarities[cacheName]["overlaps"].update(localDict[cacheName]["overlaps"])
            for otherCache, similarity in localDict[cacheName]["similarities"].items():
                cachedSimilarities[cacheName]["similarities"][otherCache] = similarity

def clusterInterlockChimera(minOvl,readIDict,maxID,minLen,minSim,contextDepths,nPairs,threadNumber):
    prevLen=0
    allReadIDict={}
    allReadIDict.update(readIDict)
    #Initializing cachedCluster
    print("Initializing cachedCluster")

    global cachedSimilarities
    cachedSimilarities=initializeCache(readIDict,"cluster",contextDepths)

    #Filling cachedCluster with similarity info
    print("Filling cachedCluster with similarity information")

    batchSize=50000

    #Create a queue for all the combinations
    combinationQueue = multiprocessing.Queue(maxsize=3*threadNumber)
    localSimilarityDictQueue = multiprocessing.Queue()
    cachedSimilarityQueue=multiprocessing.Queue()

    #Start the cache filling processes up before filling up the queue
    allProcesses=startCacheFillers(combinationQueue, localSimilarityDictQueue, cachedSimilarityQueue, cachedSimilarities, minOvl, minSim, threadNumber)
    allCacheFillers=allProcesses[0:-1]
    cacheMergerProcess=allProcesses[-1]


    #Fill up the queue with combinations
    currentBatch=[]
    i=0
    for combination in combinations(cachedSimilarities.keys(),2):
        currentBatch.append(combination)
        if len(currentBatch)>batchSize:
            combinationQueue.put(currentBatch)
            i+=1
            currentBatch=[]

    #Fill the queue with end signals
    for thread in range(0, threadNumber):
        combinationQueue.put("end")

    #Be sure to end all processes
    for cacheFiller in allCacheFillers:
        cacheFiller.join()

    localSimilarityDictQueue.put("end")

    cachedSimilarities=cachedSimilarityQueue.get()

    cacheMergerProcess.join()

    print("Preparing initial similarity index")
    bannedClusterNames=set()
    similarityIndex=generateSimilarityIndex(cachedSimilarities)
    print("Starting clustering loop ("+str(len(cachedSimilarities.keys()))+" sequences)")
    i=0
    #Parallelizing nPhase is very difficult because we cannot create processes with shared access to the cachedSimilarities dict in python.
    #Python can't yet take advantage of Copy On Write optimizations because accessing an object increases its refcount, which counts as a "write".
    #Therefore python is, in practice, Copy On Read, and therefore parallelization requires sending the updated cachedSimilarities dict to processes.
    #Sending such a large dict is prohibitively slow and ruins any improvements we might get from parallelization.
    #Future versions of python may solve this issue, but it seems doubtful. Maybe there's a different data structure I could be using.
    #If you're reading this and have any ideas, send me an email at omaroakheart@gmail.com
    #For more info on python's Copy On Read behavior look up how Instagram circumvented this problem by compiling their own custom version of CPython.
    #
    #Given the impossibility of true parallelization, I decided to at least include an optional heuristic:
    #Allow nPhase to merge more than one cluster per iterative loop.
    #Each iterative loop: Merge together the n (nPairs) non-overlapping best matches.
    #It has a low negative effect on accuracy and low positive effect on runtime in my limited testing.
    #Testing the effect of this parameters on accuracy, contiguity and runtime on large genomes is needed.
    #If you're interested in this, check out my other repo: https://github.com/OmarOakheart/Phasing-Toolkit
    while len(cachedSimilarities.keys())!=prevLen:
        prevLen=len(cachedSimilarities.keys())
        pairs,indexes=nBestPairs(cachedSimilarities,similarityIndex,minSim,maxID,bannedClusterNames,nPairs)
        for indx in reversed(indexes):
            similarityIndex.pop(indx)
        for pair in pairs:
            bannedClusterNames.update(pair)
            if pair!=[]:
                first=pair[0]
                second=pair[1]
                newName="mergedCluster_"+str(i)
                newNames=cachedSimilarities[first]["names"]|cachedSimilarities[second]["names"]
                newCluster=cachedSimilarities[first]["cluster"]|cachedSimilarities[second]["cluster"]
                commonPositions=set([x.split("=")[0] for x in cachedSimilarities[first]["consensus"]])&set([x.split("=")[0] for x in cachedSimilarities[second]["consensus"]])
                newClusterID=identity(cachedSimilarities[first]["clusterID"],cachedSimilarities[second]["clusterID"],commonPositions)
                newOverlaps=cachedSimilarities[first]["overlaps"]|cachedSimilarities[second]["overlaps"]
                newOverlaps.remove(first)
                newOverlaps.remove(second)
                for overlap in newOverlaps:
                    cachedSimilarities[overlap]["overlaps"].add(newName)
                newConsensus=consensus(newClusterID)
                newSNPPositions=set([SNP.split("=")[0] for SNP in newConsensus])
                del cachedSimilarities[first]
                del cachedSimilarities[second]

                newCache = {"names": newNames, "cluster": newCluster, "clusterID": newClusterID, "consensus": newConsensus, "positions": newSNPPositions, "similarities": {}, "overlaps": newOverlaps}
                for cache in newCache["overlaps"]:
                    cacheA = cachedSimilarities[cache]
                    cacheB = newCache
                    if first in cacheA["similarities"]:
                        del cacheA["similarities"][first]
                    if second in cacheA["similarities"]:
                        del cacheA["similarities"][second]
                    if first in cacheA["overlaps"]:
                        cacheA["overlaps"].remove(first)
                    if second in cacheA["overlaps"]:
                        cacheA["overlaps"].remove(second)
                    newCommonPos = cacheA["positions"] & cacheB["positions"]
                    if len(newCommonPos) > minOvl:
                        similarity = getSimilarity(cacheA["consensus"], cacheB["consensus"], newCommonPos, minOvl, minSim)
                    else:
                        similarity = 0
                    if similarity >= minSim:
                        if len(cacheA["consensus"]) < len(cacheB["consensus"]):
                            cacheA["similarities"][newName] = similarity
                            similarityIndex.add([similarity, len(newCommonPos), cache, newName])
                        else:
                            cacheB["similarities"][cache] = similarity
                            similarityIndex.add([similarity, len(newCommonPos), newName, cache])

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
    change=True
    while change==True:
        change=False
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

def nPhase(longReadSNPAssignments,strainName,contextDepthsFilePath,outFolder,mainFolder,referenceFilePath,minSim,minOvl,minLen,maxID,nPairs,threadNumber):
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

    baseClusters,fullClusters=clusterInterlockChimera(minOvl,readIDict,maxID,minLen,minSim,contextDepths,nPairs,threadNumber)
    baseClusters=baseClusters.values()

    baseClusterEdges=getBaseClusterEdges(100,baseClusters)

    mergedClusterEdges=mergeClusterEdges(baseClusterEdges)

    splitReadPositions=getSplitReadPositions(chimericReadIDict)
    splitReadProfiles=getSplitEdgeOverlapProfiles(splitReadPositions,mergedClusterEdges)

    usefulSplitReadProfiles=set()
    for profile in splitReadProfiles:
        if len(profile)>=2:
            usefulSplitReadProfiles.add(profile)

    ##############################################################################################################
    #TODO: Test effect of this additional cluster edge merge on accuracy?

    #mergedSplitReadProfiles=mergeClusterEdges(usefulSplitReadProfiles)
    #newChimericDict = keepUsefulSplitReadsChr(mergedSplitReadProfiles, chimericReadIDict, mergedClusterEdges)
    ##############################################################################################################

    newChimericDict=keepUsefulSplitReadsChr(usefulSplitReadProfiles,chimericReadIDict,mergedClusterEdges)

    cleanAllReadIDict={}
    cleanAllReadIDict.update(readIDict)
    cleanAllReadIDict.update(newChimericDict)

    cleanAllClusters,cleanFullClusters=clusterInterlockChimera(minOvl,cleanAllReadIDict,maxID,minLen,minSim,contextDepths,nPairs,threadNumber)

    #############
    #OUTPUT CODE#
    #############

    visDataTextFull=nPhaseFunctions.giveMeFullData(cleanAllClusters)
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
        if len(chrms)>0:
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
    minCov=nPhaseFunctions.generateCoverage(clusterReadFilePath,readFilePath,outPath,windowSize)

    outPath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_discordanceVis.tsv")
    nPhaseFunctions.generateDiscordance(clusterReadFilePath,readFilePath,outPath)

    outPath=os.path.join(outFolder,strainName+"_"+str(minOvl)+"_"+str(minSim)+"_"+str(maxID)+"_"+str(minLen)+"_minCov="+str(minCov)+"_filterVis.tsv")
    nPhaseFunctions.filterCoverage(clusterReadFilePath,readFilePath,outPath,minCov)

    return "Phasing over"



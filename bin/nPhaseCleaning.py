import itertools
import sortedcontainers
import random
import os
import glob
import bin.nPhasePipelineFunctions as nPhaseFunctions


def nPhaseCleaning(nPhaseResultFolder, longReadFilePath, strainPrefix, percentKept=98, maxDiscordance=0, deduplicate=True, FFBool=False):
    percentKept = 1 - (percentKept / 100)

    readClusterPath = glob.glob(os.path.join(nPhaseResultFolder, "Phased", "*clusterReadNames.tsv"))[0]
    knownHetFile = glob.glob(os.path.join(nPhaseResultFolder, "VariantCalls", "shortReads", "*.hetSNPs.vcf"))[0]
    readSNPsPath = glob.glob(os.path.join(nPhaseResultFolder, "VariantCalls", "longReads", "*.hetPositions.SNPxLongReads.validated.tsv"))[0]
    outPath = os.path.join(nPhaseResultFolder, "Phased", "Cleaned", strainPrefix)
    os.makedirs(outPath, exist_ok=True)
    os.makedirs(os.path.join(outPath, "cleanFastQ"), exist_ok=True)
    os.makedirs(os.path.join(outPath, "Plots"), exist_ok=True)

    phasedSNPFilePath = os.path.join(outPath, strainPrefix + "_cleaned.variants.tsv")
    phasedClusterPath = os.path.join(outPath, strainPrefix + "_cleaned.clusterReadNames.tsv")

    clusterReads = readClusterReads(readClusterPath)
    readSNPs = readReadSNPs(readSNPsPath)

    clusters = rebuildClusters(clusterReads, readSNPs)
    print("Number of clusters before cleaning:", len(clusters))

    if FFBool:
        if percentKept < 1.0:
            stitchedClusters = filterLowCoverageClusters(clusters, percentKept)
            stitchedClusters = iterativeStitching(stitchedClusters, maxDiscordance)
        else:
            stitchedClusters = iterativeStitching(clusters, maxDiscordance)
    else:
        stitchedClusters = iterativeStitching(clusters, maxDiscordance)
        if percentKept < 1.0:
            stitchedClusters = filterLowCoverageClusters(stitchedClusters, percentKept)

    if deduplicate:
        stitchedClusters = identifyGaps(stitchedClusters, readSNPs, 50)

    makeVariantsFile(stitchedClusters, phasedSNPFilePath)
    makeClusterReadFile(stitchedClusters, phasedClusterPath)

    simplerStitchedClusters = {}
    for clusterName, clusterData in stitchedClusters.items():
        simplerStitchedClusters[clusterName] = clusterData["SNPs"]

    dataText = nPhaseFunctions.giveMeFullData(simplerStitchedClusters)
    dataVisPath = os.path.join(outPath, strainPrefix + "_phasedDataFull.tsv")
    dataVisFile = open(dataVisPath, "w")
    dataVisFile.write(dataText)
    dataVisFile.close()

    simpleOutPath = os.path.join(outPath, strainPrefix + "_phasedDataSimple.tsv")
    nPhaseFunctions.simplifyDataVis(dataVisPath, simpleOutPath, 1000)
    datavisPath = os.path.join(outPath, "Plots", strainPrefix + "_")
    nPhaseFunctions.generatePhasingVis(simpleOutPath, datavisPath)

    simpleOutPath = os.path.join(outPath, strainPrefix + "_covVis.tsv")
    windowSize = 10000
    minCov = nPhaseFunctions.generateCoverage(phasedClusterPath, readSNPsPath, simpleOutPath, windowSize)

    nPhaseFunctions.generateCoverageVis(simpleOutPath, datavisPath)

    # Discordance
    simpleOutPath = os.path.join(outPath, strainPrefix + "_discordanceVis.tsv")
    nPhaseFunctions.generateDiscordance(phasedClusterPath, readSNPsPath, simpleOutPath)

    nPhaseFunctions.generateDiscordanceVis(simpleOutPath, datavisPath)

    fastQOut = os.path.join(outPath, "cleanFastQ", strainPrefix)
    nPhaseFunctions.generateLongReadFastQFiles(phasedClusterPath, longReadFilePath, fastQOut)

    print("Cleaning done! Your cleaning results are in", outPath)


def readClusterReads(filePath):
    clusterReadFile = open(filePath, "r")

    readClusterDict = {}

    for line in clusterReadFile:
        line = line.strip("\n")
        line = line.split("\t")
        readClusterDict.setdefault(line[0], set())
        readClusterDict[line[0]].add(line[1])

    return readClusterDict


def readReadSNPs(filePath):
    readSNPFile = open(filePath, "r")

    readSNPDict = {}

    for line in readSNPFile:
        line = line.strip("\n")
        line = line.split("\t")
        readID = line[0]
        for SNP in line[1:]:
            readSNPDict.setdefault(readID, set())
            readSNPDict[readID].add(SNP)

    return readSNPDict


def rebuildClusters(readClusterDict, readSNPDict):
    rebuiltClusters = {}

    for cluster, reads in readClusterDict.items():
        rebuiltClusters.setdefault(cluster, {"SNPs": {}, "Reads": set()})
        for read in reads:
            for SNP in readSNPDict[read]:
                rebuiltClusters[cluster]["SNPs"].setdefault(SNP, 0)
                rebuiltClusters[cluster]["SNPs"][SNP] += 1
                rebuiltClusters[cluster]["Reads"].add(read)

    return rebuiltClusters


def mergeClusters(clusterA, clusterB):
    mergedCluster = {"SNPs": {}, "Reads": set()}

    for cluster in [clusterA, clusterB]:
        for SNP, coverage in cluster["SNPs"].items():
            mergedCluster["SNPs"].setdefault(SNP, 0)
            mergedCluster["SNPs"][SNP] += coverage
        mergedCluster["Reads"].update(cluster["Reads"])

    return mergedCluster


def calculateDiscordance(cluster):
    discordanceScore = 0  # As a %age in the end, but weighted by coverage for each pos
    totalScore = 0

    discordanceDict = {}

    for SNP, coverage in cluster.items():
        pos = SNP.split("=")[0]
        base = SNP.split("=")[1]
        discordanceDict.setdefault(pos, {})
        discordanceDict[pos].setdefault(base, 0)
        discordanceDict[pos][base] += coverage

    for position, bases in discordanceDict.items():
        total = 0
        highest = 0
        for base, coverage in bases.items():
            total += coverage
            if coverage > highest:
                highest = coverage
        baseDiscordance = total - highest
        discordanceScore += baseDiscordance
        totalScore += total

    return discordanceScore / totalScore


def getMeanDiscordance(clusters):
    totalClusters = len(clusters)
    totalDiscordance = 0
    for clusterName in clusters:
        totalDiscordance += calculateDiscordance(clusters[clusterName]["SNPs"])
    return totalDiscordance / totalClusters


def identifyOverlap(clusterA, clusterB):
    clusterAPositions = set([SNP.split("=")[0] for SNP in clusterA.keys()])
    clusterBPositions = set([SNP.split("=")[0] for SNP in clusterB.keys()])
    return clusterAPositions & clusterBPositions


def iterativeStitching(clusters, maxDiscordance):
    readClusters = {}
    minOverlap = 10
    clusterNumber = 1

    meanDiscordance = getMeanDiscordance(clusters)
    print("Mean discordance is", meanDiscordance)

    if maxDiscordance == 0:
        maxDiscordance = meanDiscordance

    print("Initializing similarity index")

    similarityIndex = sortedcontainers.SortedList()

    for clusterA, clusterB in itertools.combinations(clusters, 2):
        overlappingRegion = identifyOverlap(clusters[clusterA]["SNPs"], clusters[clusterB]["SNPs"])
        overlapLength = len(overlappingRegion)
        if overlapLength > minOverlap:
            mergedCluster = mergeClusters(clusters[clusterA], clusters[clusterB])
            score = calculateDiscordance(mergedCluster["SNPs"])
            if score < maxDiscordance:
                similarityIndex.add([score, -overlapLength, (clusterA, clusterB)])

    print("Similarity index ready")

    print("#Number of clusters\tMerged discordance level\tInitial mean discordance\tDistance from mean\tMerged overlap length")

    while len(similarityIndex) > 0:
        highestOverlapLength = similarityIndex[0][1]
        lowestScore = similarityIndex[0][0]
        lowestPair = similarityIndex[0][2]
        if lowestScore < maxDiscordance:
            clusters["stitched_" + str(clusterNumber)] = mergeClusters(clusters[lowestPair[0]], clusters[lowestPair[1]])
            del clusters[lowestPair[0]]
            del clusters[lowestPair[1]]
            getRidOf = []
            for pair in similarityIndex:
                if lowestPair[0] in pair[2] or lowestPair[1] in pair[2]:
                    getRidOf.append(pair)
            for pair in getRidOf:
                similarityIndex.discard(pair)
            for cluster in clusters.keys():
                if cluster != "stitched_" + str(clusterNumber):
                    overlappingRegion = identifyOverlap(clusters[cluster]["SNPs"],
                                                        clusters["stitched_" + str(clusterNumber)]["SNPs"])
                    overlapLength = len(overlappingRegion)
                    if overlapLength > minOverlap:
                        newMergedCluster = mergeClusters(clusters[cluster], clusters["stitched_" + str(clusterNumber)])
                        newScore = calculateDiscordance(newMergedCluster["SNPs"])
                        if newScore < maxDiscordance:
                            similarityIndex.add([newScore, -overlapLength, (cluster, "stitched_" + str(clusterNumber))])  ##
            clusterNumber += 1
        print(len(clusters), round(lowestScore * 100, 2), round(meanDiscordance * 100, 2),
              round(meanDiscordance * 100 - lowestScore * 100, 2), -highestOverlapLength)

    print("Merging step over")

    return clusters


def filterLowCoverageClusters(clusters, percentKept):
    # calculate coverage per cluster
    clusterCoverageList = []
    totalCoverage = 0
    for cluster in clusters.keys():
        clusterCoverage = 0
        for SNP, coverage in clusters[cluster]["SNPs"].items():  # Maybe only count reads that support the consensus?
            clusterCoverage += coverage
        # calculate total coverage
        totalCoverage += clusterCoverage
        clusterCoverageList.append((clusterCoverage, cluster))
    # sort clusters by coverage
    clusterCoverageList.sort(reverse=True)
    # take most covered clusters until you reach 95% of total coverage
    keepN = percentKept * totalCoverage
    filteredClusters = {}
    filteredCoverage = 0
    for pair in clusterCoverageList:
        cov = pair[0]
        clust = pair[1]
        if filteredCoverage < keepN:
            filteredClusters[clust] = clusters[clust]
            filteredCoverage += cov
        else:
            return filteredClusters


def identifyGaps(filteredClusters, readSNPs, minGapLen):
    redistributedClusters = {}
    # Start by identifying all positions and sort them by chromosome
    genomeCoverage = {}
    for cluster in filteredClusters.keys():
        theseSNPs = filteredClusters[cluster]["SNPs"].keys()
        uniquePositions = set([SNP.split("=")[0] for SNP in theseSNPs])
        for pos in uniquePositions:
            chr = pos.split(":")[0]
            pos = int(pos.split(":")[1])
            genomeCoverage.setdefault(chr, {})
            genomeCoverage[chr].setdefault(pos, 0)
            genomeCoverage[chr][pos] += 1
    # What is a gap?
    for chromosome in genomeCoverage.keys():
        redistributedClusters.setdefault(chromosome, {})
        chromosomePositions = []
        for pos, cov in genomeCoverage[chromosome].items():
            chromosomePositions.append((pos, cov))
        chromosomePositions.sort()
        # Identify a gap in here
        ploidy = round(sum([cov for pos, cov in chromosomePositions]) / len(chromosomePositions))
        gapPositionDict = {}
        i = 1
        for position, cov in chromosomePositions:
            if cov < ploidy:
                gapPositionDict.setdefault("gap_" + str(i), set())
                gapPositionDict["gap_" + str(i)].add(chromosome + ":" + str(position))
            else:
                i += 1
        tinyGapNames = set()
        for gap in gapPositionDict.keys():
            gapLength = len(gapPositionDict[gap])
            if gapLength < minGapLen:
                tinyGapNames.add(gap)
        for tinyGap in tinyGapNames:
            del gapPositionDict[tinyGap]

        for gap, allGaps in gapPositionDict.items():
            allGaps = [int(SNP.split(":")[1]) for SNP in allGaps]
            allGaps.sort()

        # Okay now we have gaps of adequate size
        # Extract the clusters in this region and their coverage levels.
        clusterPositionDict = {}
        for cluster in filteredClusters.keys():
            clusterPositions = []
            for SNP, coverage in filteredClusters[cluster]["SNPs"].items():
                clusterPositions.append(SNP.split("=")[0])
            clusterPositionDict[cluster] = set(clusterPositions)
        # Now identify, for each gap, every single cluster that overlaps with it.
        gapOverlapDict = {}
        for gap, gapPositions in gapPositionDict.items():
            gapOverlapDict.setdefault(gap, [])
            mergeCandidates = []
            for cluster, clusterPositions in clusterPositionDict.items():
                overlap = gapPositions & clusterPositions
                if overlap != set():
                    gapOverlapDict[gap].append((len(overlap), cluster))
                if overlap == set() and list(clusterPositions)[0].split(":")[
                    0] == chromosome:  # If the cluster's in the right chromosome, doesn't overlap the gap
                    mergeCandidates.append((len(overlap), cluster))
            gapOverlapDict[gap].sort(reverse=True)
            mergeCandidates.sort(reverse=True)
            # Now, for each cluster that overlaps a gap, let's try redistributing reads.
            # 1, identify which of the donor clusters has the highest coverage value for the region
            highestCoverage = 0
            highestCandidate = None
            for overlapLen, donorCandidate in gapOverlapDict[gap]:
                if overlapLen > 0:
                    coverage = getRegionCoverage(filteredClusters[donorCandidate]["SNPs"], gapPositions)
                    if coverage > highestCoverage:
                        highestCoverage = coverage
                        highestCandidate = donorCandidate
            # 1.5, identify the reads of that higher coverage cluster
            includedReads, remainingReads = getRegionReads(filteredClusters[highestCandidate], gapPositions, readSNPs)
            # 2, split the reads in that cluster and region into at least two different clusters (and minimize discordance with 2 groups)
            if len(includedReads) > 20 and len(gapOverlapDict[gap]) > 0 and len(mergeCandidates) > 0:
                del filteredClusters[highestCandidate]
                clusterA, clusterB, readsA, readsB = minimizeRegionDiscordance(includedReads, readSNPs)

                # 3, figure out the lowest discordance merge you can do (with the other clusters of the chromosome
                # Make sure you don't double count reads when you overlap several gaps.

                # Identify which of clusterA and clusterB merges best back into the original cluster
                levelA, mergedA = checkDiscordance(readsA | remainingReads, readSNPs)
                levelB, mergedB = checkDiscordance(readsB | remainingReads, readSNPs)

                if levelA < levelB:
                    filteredClusters[highestCandidate] = mergedA
                    donorCluster = makeCluster(readsB, readSNPs)  # Or just clusterB
                else:
                    filteredClusters[highestCandidate] = mergedB
                    donorCluster = makeCluster(readsA, readSNPs)  # Or just clusterA

                # Replace the donor region with the reads that were chosen

                lowestDiscordance = 1000
                lowestCandidate = None
                lowestCluster = None

                for candidateLength, candidate in mergeCandidates:
                    # Test every combination of candidate and remaining cluster

                    # Calculate discordance
                    readClusterToTest = donorCluster["Reads"] | filteredClusters[candidate]["Reads"]
                    testDiscordance, mergedDonor = checkDiscordance(readClusterToTest, readSNPs)
                    if testDiscordance < lowestDiscordance:
                        lowestDiscordance = testDiscordance
                        lowestCandidate = candidate
                        lowestCluster = mergedDonor

                # Keep the one with the lowest discordance
                # Merge the one with the lowest discordance
                del filteredClusters[lowestCandidate]
                filteredClusters[lowestCandidate] = lowestCluster

            else:
                pass
    return filteredClusters


def checkDiscordance(readsFound, readSNPs):
    SNPCluster = {"Reads": set(), "SNPs": {}}
    for read in readsFound:
        for SNP in readSNPs[read]:
            SNPCluster["SNPs"].setdefault(SNP, 0)
            SNPCluster["SNPs"][SNP] += 1
        SNPCluster["Reads"].add(read)
    return calculateDiscordance(SNPCluster["SNPs"]), SNPCluster


def getRegionReads(cluster, allowedPositions, readSNPs):
    includedReads = set()
    remainingReads = set()
    for read in cluster["Reads"]:
        allSNPs = list(readSNPs[read])
        allPositions = set([SNP.split("=")[0] for SNP in allSNPs])
        if len(allPositions & allowedPositions) > 5:
            includedReads.add(read)
        else:
            remainingReads.add(read)
    return includedReads, remainingReads


def makeCluster(reads, readSNPs):
    SNPCluster = {"Reads": set(), "SNPs": {}}
    for read in reads:
        for SNP in readSNPs[read]:
            SNPCluster["SNPs"].setdefault(SNP, 0)
            SNPCluster["SNPs"][SNP] += 1
        SNPCluster["Reads"].add(read)
    return SNPCluster


def minimizeRegionDiscordance(includedReads, readSNPs):
    startingCluster = {}
    consensusSNPs = set()

    # Identify one read that feels like it belongs to cluster #2
    for read in includedReads:
        for SNP in readSNPs[read]:
            pos = SNP.split("=")[0]
            base = SNP.split("=")[1]
            startingCluster.setdefault(pos, {})
            startingCluster[pos].setdefault(base, 0)
            startingCluster[pos][base] += 1

    for position, positionData in startingCluster.items():
        maxCov = 0
        maxBase = None
        for possibleBase, cov in positionData.items():
            if cov > maxCov:
                maxCov = cov
                maxBase = possibleBase
        consensusSNPs.add(position + "=" + maxBase)

    # Identify one read that feels like it belongs to cluster #2
    # Find two reads that disagree the most (even if they disagree very little)
    maxNotInCommon = 0
    maxPair = (None, None)

    for readA, readB in itertools.combinations(includedReads, 2):
        SNPsA = readSNPs[readA]
        SNPsB = readSNPs[readB]
        posA = set([SNP.split("=")[0] for SNP in SNPsA])
        posB = set([SNP.split("=")[0] for SNP in SNPsB])
        commonPos = posA & posB
        commonSNPs = SNPsA & SNPsB
        SNPsNotInCommon = len(commonPos) - len(commonSNPs)
        if SNPsNotInCommon >= maxNotInCommon:
            maxNotInCommon = SNPsNotInCommon
            maxPair = (readA, readB)
    # Grow both clusters
    clusterA = {}
    clusterB = {}
    alreadyUsed = {maxPair[0], maxPair[1]}
    for SNP in readSNPs[maxPair[0]]:
        pos = SNP.split("=")[0]
        base = SNP.split("=")[1]
        clusterA.setdefault(pos, {})
        clusterA[pos].setdefault(base, 0)
        clusterA[pos][base] += 1

    readsA = set()
    readsB = set()

    while len(alreadyUsed) != len(includedReads):
        clusterAConsensus = getConsensus(clusterA)
        clusterBConsensus = getConsensus(clusterB)
        bestA = None
        bestAScore = 0
        bestB = None
        bestBScore = 0
        for read in includedReads:
            if read not in alreadyUsed:
                candidateSNPs = readSNPs[read]
                aScore = len(candidateSNPs & clusterAConsensus)
                if aScore >= bestAScore:
                    bestAScore = aScore
                    bestA = read
                bScore = len(candidateSNPs & clusterBConsensus)
                if bScore >= bestBScore:
                    bestBScore = bScore
                    bestB = read
        if bestB == bestA:
            x = random.choice([0, 1])
            if x == 0:
                for SNP in readSNPs[bestA]:
                    readsA.add(bestA)
                    pos = SNP.split("=")[0]
                    base = SNP.split("=")[1]
                    clusterA.setdefault(pos, {})
                    clusterA[pos].setdefault(base, 0)
                    clusterA[pos][base] += 1
                alreadyUsed.add(bestA)
            elif x == 1:
                for SNP in readSNPs[bestB]:
                    readsB.add(bestB)
                    pos = SNP.split("=")[0]
                    base = SNP.split("=")[1]
                    clusterB.setdefault(pos, {})
                    clusterB[pos].setdefault(base, 0)
                    clusterB[pos][base] += 1
                alreadyUsed.add(bestB)
        else:
            if bestA is not None:
                readsA.add(bestA)
                for SNP in readSNPs[bestA]:
                    pos = SNP.split("=")[0]
                    base = SNP.split("=")[1]
                    clusterA.setdefault(pos, {})
                    clusterA[pos].setdefault(base, 0)
                    clusterA[pos][base] += 1
                alreadyUsed.add(bestA)
            if bestB is not None:
                readsB.add(bestB)
                for SNP in readSNPs[bestB]:
                    pos = SNP.split("=")[0]
                    base = SNP.split("=")[1]
                    clusterB.setdefault(pos, {})
                    clusterB[pos].setdefault(base, 0)
                    clusterB[pos][base] += 1
                alreadyUsed.add(bestB)

    # Do that iteratively until you have two different clusters and nothing moves anymore
    return clusterA, clusterB, readsA, readsB


def getConsensus(cluster):
    consensusSequence = set()
    for SNP, bases in cluster.items():
        maxBase = None
        maxCov = 0
        for base, cov in bases.items():
            if cov > maxCov:
                maxCov = cov
                maxBase = base
        consensusSequence.add(maxBase)
    return consensusSequence


def getRegionCoverage(cluster, allowedPositions):
    totalCoverage = 0  # As a %age in the end, but weighted by coverage for each pos

    for SNP, coverage in cluster.items():
        pos = SNP.split("=")[0]
        if pos in allowedPositions:
            totalCoverage += coverage

    return totalCoverage / len(allowedPositions)


def getRegionDiscordance(cluster, allowedPositions):
    discordanceScore = 0  # As a %age in the end, but weighted by coverage for each pos
    totalScore = 0

    discordanceDict = {}

    for SNP, coverage in cluster.items():
        pos = SNP.split("=")[0]
        if pos in allowedPositions:
            base = SNP.split("=")[1]
            discordanceDict.setdefault(pos, {})
            discordanceDict[pos].setdefault(base, 0)
            discordanceDict[pos][base] += coverage

    for position, bases in discordanceDict.items():
        total = 0
        highest = 0
        for base, coverage in bases.items():
            total += coverage
            if coverage > highest:
                highest = coverage
        baseDiscordance = total - highest
        discordanceScore += baseDiscordance
        totalScore += total

    if totalScore != 0:
        return discordanceScore / totalScore
    else:
        return 1000000


# Then basically take the regions that are too highly covered and redistribute their reads to the best cluster
# Either it stays in its own cluster, or it goes in another, just pick the best option each time
# That means identify the read that fucks the most with the system, then remove it, then see where it fits best
# Once a read has been moved, it's "safe" and even if placed back in the cluster gets ignored
# Then in the end either you get rid of a cluster or you make it more likely to be stitched back in.

def makeVariantsFile(stitchedClusters, filePath):
    variantFileText = ""
    for cluster, clusterData in stitchedClusters.items():
        posDict = {}
        for SNP, coverage in clusterData["SNPs"].items():
            chromosome = SNP.split(":")[0]
            position = SNP.split(":")[1].split("=")[0]
            base = SNP.split("=")[1]
            posDict.setdefault(chromosome + ":" + position, [0, None])
            if coverage > posDict[chromosome + ":" + position][0]:
                posDict[chromosome + ":" + position][0] = coverage
                posDict[chromosome + ":" + position][1] = base
        for position, data in posDict.items():
            chromosome = position.split(":")[0]
            pos = position.split(":")[1]
            base = data[1]
            variantFileText += "\t".join([cluster, chromosome, pos, base]) + "\n"

    variantFile = open(filePath, "w")
    variantFile.write(variantFileText)
    variantFile.close()

def makeClusterReadFile(stitchedClusters, filePath):
    variantFileText = ""
    for cluster, clusterData in stitchedClusters.items():
        for read in clusterData["Reads"]:
            variantFileText += cluster + "\t" + read + "\n"
    variantFile = open(filePath, "w")
    variantFile.write(variantFileText)
    variantFile.close()

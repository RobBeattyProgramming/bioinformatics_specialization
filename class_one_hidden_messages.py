import math

def PatternCount(text, pattern):
    count = 0
    kLength = len(pattern)
    
    for i in range(0, len(text) - kLength + 1):
        if text[i:i + kLength] == pattern:
            count += 1
    return count

def FrequentWords(text, kmer):
    kmerCount = {}
    maxKmerDict = {}
 
    for i in range(0, len(text) - kmer + 1):
        kmerCount[text[i:i + kmer]] = PatternCount(text, text[i:i + kmer])
    
    maxKmerValue = max(kmerCount.values())

    for i in kmerCount:
        if kmerCount[i] == maxKmerValue:
            maxKmerDict[i] = maxKmerValue

    return maxKmerDict

def ReverseComplement(dnaString):
    dnaString = dnaString[::-1]
    reverseComplement = ""

    for i in dnaString:
        if i == "A":
            reverseComplement = reverseComplement + "T"
        if i == "T":
            reverseComplement = reverseComplement + "A"
        if i == "C":
            reverseComplement = reverseComplement + "G"
        if i == "G":
            reverseComplement = reverseComplement + "C"
    
    return reverseComplement

def PatternMatching(pattern, genome):
    complementPattern = ReverseComplement(pattern)
    patternPoints = []
    count = 0

    for i in range(0, len(genome) - len(pattern) + 1):
        if genome[i:i + len(pattern)] == pattern:
            patternPoints.append(count)
        elif genome[i:i + len(pattern)] == complementPattern:
            patternPoints.append(count)
        count += 1

    return patternPoints

def FindClumps(kmerSize, clumpLength, repetitionAmount, genome):
    kmerDict = {}
    qualifiedKmerDict = {}
    finalDict = {} 
    for i in range(0, len(genome) - clumpLength + 1):
        currentKmer = genome[i:i + kmerSize]
        if currentKmer in kmerDict:
            kmerDict[currentKmer] = kmerDict[currentKmer] + 1
        else:
            kmerDict[currentKmer] = 1
    for key in kmerDict:
        if kmerDict[key] >= repetitionAmount:
            qualifiedKmerDict[key] = kmerDict[key]
    for key in qualifiedKmerDict:
        currentKmerRepetition = PatternMatching(key, genome)
        currentKmerRepetition.reverse()
        clumpRange = 1 
        for i in currentKmerRepetition:
            for j in range(1, len(currentKmerRepetition) - 1):
                if i - currentKmerRepetition[j] <= clumpLength:
                    clumpRange = clumpRange + 1
                else:
                    break
                
        if clumpRange >= repetitionAmount:
            finalDict[key] = currentKmerRepetition
    return finalDict

# Module 2

def Skew(genome):
    skew = [0]
    currentSkew = 0

    for i in genome:
        if i == "C":
            currentSkew = currentSkew - 1
        elif i == "G":
            currentSkew = currentSkew + 1
        skew.append(currentSkew)

    return skew

def MinimumSkew(genome):
    skewGraph = Skew(genome)
    minimumValue = min(skewGraph)

    minimumLocations = []
    currentLocation = 0
    for i in skewGraph:
        if i == minimumValue:
            minimumLocations.append(currentLocation)
        currentLocation = currentLocation + 1
    
    return minimumLocations

def HammingDistance(p,q):
    qDistance = 0
    hammingDistanceAmount = 0
    for i in p:
        if i == q[qDistance]:
            pass
        else:
            hammingDistanceAmount = hammingDistanceAmount + 1

        qDistance = qDistance + 1

    return hammingDistanceAmount

def ApproximatePatternMatching(text, genome, hd):
    startingPosition = 0
    startingPositionLocations = []
    for i in range(len(genome) - len(text) + 1):
        if genome[i:i+len(text)] == text:
            startingPositionLocations.append(startingPosition)
            pass
        elif HammingDistance(text, genome[i:i+len(text)]) <= hd:
            startingPositionLocations.append(startingPosition)
            
        startingPosition = startingPosition + 1

    return startingPositionLocations

def PatternCountTwo(genome, text, hd):
    count = 0
    for i in range(len(genome) - len(text) + 1):
        if genome[i:i+len(text)] == text:
            count = count + 1
            pass
        elif HammingDistance(text, genome[i:i+len(text)]) <= hd:
            count = count + 1
    return count

def ImmediateNeighbors(pattern):
    neighborhood = []

    for i in range(0, len(pattern)):
        nucleotides = ["A", "C", "G", "T"]
        patternList = list(pattern)
        currentPatternNucleotide = patternList[i]
        nucleotides.remove(currentPatternNucleotide)
         
        for j in nucleotides:
            patternList[i] = j
            currentNeighbor = "".join(patternList)
            neighborhood.append(currentNeighbor)

            patternList = list(pattern)

    neighborhood.append(pattern)

    return neighborhood

def Neighbors(pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]
    
    dOneNeighborhood = ImmediateNeighbors(pattern)

    if d == 1:
        return dOneNeighborhood
    else:
        neighborhood = []
        i = 1
        while i != d:
            for j in dOneNeighborhood:
                currentOneNieghborhood = ImmediateNeighbors(j)
                for z in currentOneNieghborhood:
                    neighborhood.append(z)
            i = i + 1

        for xx in neighborhood:
            hammingDistance = HammingDistance(pattern, xx)
            if hammingDistance > d:
                neighborhood.remove(xx)

        return list(set(neighborhood))

def FrequentWordsWithMismatches(text, k, d):
    patternCount = {}
    maxKmers = []
    
    for i in range(0, len(text) - k + 1):
        currentKmer = text[i:i + k]
        currentKmerNeighborhood = Neighbors(currentKmer, d)


        for l in currentKmerNeighborhood:
            if l not in patternCount:
                patternCount[l] = 1
            else:
                count = patternCount[l]
                updatedCount = count + 1
                patternCount[l] = updatedCount
    
    maxValue = max(patternCount.values())
    for x in patternCount:
        if patternCount[x] == maxValue:
            maxKmers.append(x)

    return maxKmers         

def FrequentWordsWithMismatchesAndReverseComplements(text, k, d):
    patternCount = {}
    maxKmers = []
    
    for i in range(0, len(text) - k + 1):
        currentKmer = text[i:i + k]
        reverseKmer = ReverseComplement(currentKmer)
        currentKmerNeighborhood = Neighbors(currentKmer, d)
        reverseNeighborhood = Neighbors(reverseKmer, d)


        for l in currentKmerNeighborhood:
            if l not in patternCount:
                patternCount[l] = 1
            else:
                count = patternCount[l]
                updatedCount = count + 1
                patternCount[l] = updatedCount

        for lr in reverseNeighborhood:
            if lr not in patternCount:
                patternCount[lr] = 1
            else:
                count = patternCount[lr]
                updatedCount = count + 1
                patternCount[lr] = updatedCount
    
    reverseGenome = ReverseComplement(text)

    for i in range(0, len(reverseGenome) - k + 1):
        currentKmer = text[i:i + k]
        reverseKmer = ReverseComplement(currentKmer)
        currentKmerNeighborhood = Neighbors(currentKmer, d)
        reverseNeighborhood = Neighbors(reverseKmer, d)


        for j in currentKmerNeighborhood:
            if j not in patternCount:
                patternCount[j] = 1
            else:
                count = patternCount[j]
                updatedCount = count + 1
                patternCount[j] = updatedCount

        for jr in reverseNeighborhood:
            if jr not in patternCount:
                patternCount[jr] = 1
            else:
                count = patternCount[jr]
                updatedCount = count + 1
                patternCount[jr] = updatedCount
    
    reverseGenome = ReverseComplement(text)

    
    maxValue = max(patternCount.values())
    for x in patternCount:
        if patternCount[x] == maxValue:
            maxKmers.append(x)

    return maxKmers

# Module 3  

def MotifEnumeration(dna, k, d):
    neighbors = {}
    dnaList = []

    dnaSegment = ""
    for i in dna:
        if i == " ":
            dnaList.append(dnaSegment)
            dnaSegment = ""
        else:
            dnaSegment = dnaSegment + i
    dnaList.append(dnaSegment)

    currentSegment = 1
    for segment in dnaList:
        for x in range(0, len(segment) - k + 1):
            kmer = segment[x: x + k]
            kmerNeighbors = Neighbors(kmer, d)

            for z in kmerNeighbors:
                if z not in neighbors:
                    neighbors[z] = [currentSegment]
                else:
                    currentKList = neighbors[z]
                    currentKList.append(currentSegment)
                    neighbors[z] = currentKList


        currentSegment = currentSegment + 1
    
    neededNumber = len(dnaList)
    confirmationList = []
    for c in range(1, neededNumber + 1):
        confirmationList.append(c)
    
    qualifiedKmers = []
    
    for potentialK in neighbors:
        finalCount = 0
        for kCount in confirmationList:
            if kCount in neighbors[potentialK]:
                finalCount = finalCount + 1
        if finalCount ==  neededNumber:
            qualifiedKmers.append(potentialK)
    return qualifiedKmers


test = "TCGGGGGTTTTT CCGGTGACTTAC ACGGGGATTTTC TTGGGGACTTTT AAGGGGACTTCC TTGGGGACTTCC TCGGGGATTCAT TCGGGGATTCCT TAGGGGAACTAC TCGGGTATAACC"

def profile(strings):
    
    individualStrings = []
    entropyList = []
    dnaSegment = ""
    for i in strings:
        if i == " ":
            individualStrings.append(dnaSegment)
            dnaSegment = ""
        else:
            dnaSegment = dnaSegment + i
    individualStrings.append(dnaSegment)

    for x in range(0,len(individualStrings[0]) - 1):

        nucleotideCount = {
            "A":0,
            "C":0,
            "G":0,
            "T":0
        }

        
        for y in individualStrings:
            if y[x] == "A":
                count = nucleotideCount["A"]
                updatedCount = count + 1
                nucleotideCount["A"] = updatedCount
            elif y[x] == "C":
                count = nucleotideCount["C"]
                updatedCount = count + 1
                nucleotideCount["C"] = updatedCount
            elif y[x] == "G":
                count = nucleotideCount["G"]
                updatedCount = count + 1
                nucleotideCount["G"] = updatedCount
            elif y[x] == "T":
                count = nucleotideCount["T"]
                updatedCount = count + 1
                nucleotideCount["T"] = updatedCount
        
        a = float(nucleotideCount["A"]) / len(individualStrings)
        c = float(nucleotideCount["C"]) / len(individualStrings)
        g = float(nucleotideCount["G"]) / len(individualStrings)
        t = float(nucleotideCount["T"]) / len(individualStrings)
        profileList = [a, c, g, t]

        aboveZero = []
        for profile in profileList:
            if profile > 0:
                aboveZero.append(profile)

        loggedList = []
        for pp in aboveZero:
            logged = math.log(pp, 2)
            logged = pp * logged
            loggedList.append(logged)
        
        summeedProfle = 0
        for ppp in loggedList:   
            summeedProfle = summeedProfle + ppp
        
        completedProfile = -1 * summeedProfle
        roundedProfile = round(completedProfile, 3)
        
        entropyList.append(roundedProfile)

    return entropyList

        


print(profile(test))
            

            

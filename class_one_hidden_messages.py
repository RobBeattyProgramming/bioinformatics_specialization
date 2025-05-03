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
    
 
"""def FindClumps(s, l, t, genome):
    clumpDict = {}

    for i in range(0, len(genome) - l + 1):
        repetitionPoints = PatternMatching(genome[i:i + s], genome[i:i+l])
        if len(repetitionPoints) >= t:
            clumpDict[genome[i:i + s]] = len(repetitionPoints)
    
    return clumpDict"""


"""def FindClumps(kmerSize, clumpLength, repetitionAmount, genome):
    #call if kmer in test, continue so less run time
    kmerRepetitionDict = {}

    #iterates over genome in length of kmer
    for i in range(0, len(genome) - clumpLength + 1):
        clumpDict = {} #use for more info 

        currentKmer = genome[i:i + kmerSize]
        #checks for repeats of kmer and skips over if it is already in dict
        if currentKmer in kmerRepetitionDict:
            continue
        else:
            #checks repetition list via PatternMatching to see if at least t of kmer are in range
            currentKmerRepetion = PatternMatching(currentKmer, genome)
            if len(currentKmerRepetion) < repetitionAmount:
                continue
            else:
                clumpRange = 0
                for i in currentKmerRepetion:
                    for j in range(1, len(currentKmerRepetion) - 1):
                        if i - currentKmerRepetion[j] <= clumpLength:
                            clumpRange = clumpRange + 1
                        else:
                            break
                
                if clumpRange >= repetitionAmount:
                    clumpDict[currentKmer] = currentKmerRepetion
      
    return clumpDict"""

"""#attempt of FindClumps above works correctly but is impossible to execute in large sizes
def FindClumps(kmerSize, clumpLength, repetitionAmount, genome):
    kmerDict = {}
    qualifiedKmerDict = {}
    finalDict = {} #clean up by just removing from qualified kmer dict
    #iterates genome in segments of length kmer and tallies repeats
    for i in range(0, len(genome) - clumpLength + 1):
        currentKmer = genome[i:i + kmerSize]
        if currentKmer in kmerDict:
            # NEEDS TO PUT LOCATIONS LIKE PATTERN MATCH 
            kmerDict[currentKmer] = kmerDict[currentKmer] + 1
        else:
            kmerDict[currentKmer] = 1
    #iterates tallied kmer dict and adds to qualified kmer dict if repeats are equal or larger to repetitionAmount
    for key in kmerDict:
        if kmerDict[key] >= repetitionAmount:
            qualifiedKmerDict[key] = kmerDict[key]
    #PatternMatches qualified kmers 
    for key in qualifiedKmerDict:
        currentKmerRepetition = PatternMatching(key, genome)
        currentKmerRepetition.reverse()
        clumpRange = 1 #changed from 0
        for i in currentKmerRepetition:
            for j in range(1, len(currentKmerRepetition) - 1):
                if i - currentKmerRepetition[j] <= clumpLength:
                    clumpRange = clumpRange + 1
                else:
                    break
                
        if clumpRange >= repetitionAmount:
            finalDict[key] = currentKmerRepetition

    
    return finalDict"""

#attempt of FindClumps above works correctly but is impossible to execute in large sizes
def FindClumps(kmerSize, clumpLength, repetitionAmount, genome):
    kmerDict = {}
    dupDict = {}

    currentLocation = 0
    #iterates genome in segments of length kmer and adds pattern locations
    for i in range(0, len(genome) - clumpLength + 1):
        currentKmer = genome[i:i + kmerSize]
        if currentKmer in kmerDict:
            kmerDict[currentKmer].append(currentLocation)
        else:
            kmerDict[currentKmer] = [currentLocation]
        currentLocation = currentLocation + 1

        dupDict = kmerDict
        for key in dupDict:
            if len(dupDict[key]) < repetitionAmount:
                del kmerDict[key]
                

    return kmerDict
x = open("ecoli.txt", "r")
toot = x.read()

testy = toot[0:100000]

her = FindClumps(9, 500, 3, testy)
print(her)





 
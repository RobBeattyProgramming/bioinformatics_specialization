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
    
 
def FindClumps(s, l, t, genome):
    clumpDict = {}

    for i in range(0, len(genome) - l + 1):
        repetitionPoints = PatternMatching(genome[i:i + s], genome[i:i+l])
        if len(repetitionPoints) >= t:
            clumpDict[genome[i:i + s]] = len(repetitionPoints)
    
    return clumpDict

her = open("vib.txt", "r")
x = her.read()

"""all that's happening past this is seeing if there's enough
    patterns to justify putting it into dictionary to return.
    how can this be reduced?? """

"""WAIT WAIT IS THE PROBLEM JUST THAT PATTERN MATCHING IS 
    CHECKING THE WHOLE GENOME AND NOT THE CLUMP???"""

print(FindClumps(9, 500, 3, x))

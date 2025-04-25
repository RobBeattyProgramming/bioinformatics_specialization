vibrioCholeraeOri = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"

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
    
    print()
    print(kmerCount)
    print()


    return maxKmerDict


print(FrequentWords(vibrioCholeraeOri, 9))




    




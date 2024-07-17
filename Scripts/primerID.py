#simple method, ignore indel divergence. Output a "fake" R1 for subsequent STAR mapping. Must have majority nucleotides at all positions. Set arbitrary score at Q40 (ASCII I)
from statistics import mode

def highQualityfastq(R1, R2, adaptUp, UMIlen, outFile, trimSeq):
    data = {}
    with open(R1, 'r') as R1, open(R2, 'r') as R2:
        line1 = R1.readline()
        line2 = R2.readline()
        numTraversed = 1
        while line1 != '':
            if numTraversed % 4 == 1:
                name1 = line1
                name2 = line2
                assert name1.split(' ')[0] == name2.split(' ')[0]
            elif numTraversed % 4 == 2:
                seq1 = line1
                seq2 = line2
            elif numTraversed %4 == 0:
                pos = seq1.find(adaptUp)
                if pos != -1:
                    barcode = seq1[pos + len(adaptUp):pos+len(adaptUp) + UMIlen]
                    currDict = data

                    for element in barcode:
                        if element not in currDict:
                            currDict[element] = {}
                        currDict = currDict[element]
                    if len(currDict) == 0:
                        currDict[element] = []
                    currDict[element] += [seq2[:-1]]
 
            numTraversed = (numTraversed + 1)
            line1 = R1.readline()
            line2 = R2.readline()
    output = {}
    with open(outFile, 'w') as outfile:
        numSeq = 1

        for value in list(NestedDictValues(data)):
            if len(value) >= 3:
                seq = ''
                #if any position does not have a mode throw everything away.
                fullSeq = True
                for position, toss in enumerate(value[0]):
                    baseObs = []
                    for seq in value:
                        baseObs += seq[position]
                    try:
                        base = mode(baseObs)
                        seq += base
                    except:
                        fullSeq = False
                        pass
                if fullSeq:
                    outfile.write('@Sequence' + str(numSeq) +'barcode' + str(len(value)) + '\n')
                    numSeq += 1
                    outfile.write(seq[trimSeq:] + '\n')
                    outfile.write('+\n')
                    score = ''
                    for character in seq[trimSeq:]:
                        score += 'I'
                    outfile.write(score + '\n')
            
                



#from https://stackoverflow.com/questions/23981553/get-all-values-from-nested-dictionaries-in-python
def NestedDictValues(d):
  for v in d.values():
    if isinstance(v, dict):
      yield from NestedDictValues(v)
    else:
      yield v
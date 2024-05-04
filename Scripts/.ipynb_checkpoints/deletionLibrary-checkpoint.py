'''methods for designing and analyzing deletion-spanning library'''
import math
import pandas as pd
import re
from string import ascii_uppercase

'''design primers, attempting to normalize meltemp.
Stagger according the the scaled staggering, based on the provided fasta. Returns error
if not all junctions would be unique for mapping'''
def makeDeletionPrimers(*, staggerList, oligoConcentration, sodiumConcentration,
        targetSequence, targetTemperature, minimumLength):
    #each interval consists of a start positon, followed  by the interval at which deletions will be made. Primers designed in inverse position.
    #assumes a 1-indexed system of base position, as most do
    sequenceFrame = []
    for interval in staggerList:
        firstBase = interval [0]
        lastBase = interval [1]
        stagger = interval[2]
        while firstBase <= lastBase:
            #start with 2 bases for meltemp calculation
            currLength = minimumLength
            sequence = reverseComplement(targetSequence[firstBase - currLength:firstBase])

            while(melTemp(sequence=sequence, sodiumConcentration=sodiumConcentration, oligoConcentration=oligoConcentration) < targetTemperature):
                currLength += 1
                sequence = reverseComplement(targetSequence[firstBase - currLength:firstBase])
            sequenceFrame += [pd.DataFrame({'annealingSequence':[sequence], 'inferredMeltemp':[melTemp(sequence=sequence, sodiumConcentration=sodiumConcentration, oligoConcentration=oligoConcentration)], 'position':[firstBase]})]
            firstBase += stagger

    sequenceFrame = pd.concat(sequenceFrame)
    return sequenceFrame



def reverseComplement(sequence):
    sequenceList = [character for character in sequence]
    revSeq = ''
    while len(sequenceList) > 0:
        character = sequenceList.pop()
        if character == 'A':
            revSeq += 'T'
        elif character == 'T':
            revSeq += 'A'
        elif character == 'G':
            revSeq += 'C'
        elif character == 'C':
            revSeq += 'G'
        else:
            revSeq += 'N'
    return revSeq


def melTemp(*, sequence, sodiumConcentration, oligoConcentration):
    #nearest-neighbor chart. Culled from Breslauer et al. 1986 PNAS. enthalpy
     #∆H and entropy ∆S provided.
    meltingTable = {'AA':{'H':9.1, 'S':24.0}, 'AT':{'H':8.6, 'S':23.9},
    'TA':{'H':6.0, 'S':16.9}, 'CA':{'H':5.8, 'S':12.9}, 'GT':{'H':6.5, 'S':17.3},
    'CT':{'H':7.8, 'S':20.8}, 'GA':{'H':5.6, 'S':13.5}, 'CG':{'H':11.9, 'S':27.8},
    'GC':{'H':11.1, 'S':26.7}, 'GG':{'H':11.0, 'S':26.6}}
    more = {}
    for key in meltingTable:
        revKey = reverseComplement(key)
        if revKey not in meltingTable:
            more[revKey] = meltingTable[key]
    meltingTable.update(more)
    runningEnthalpy = 0
    runningEntropy = 0
    for position, character in enumerate(sequence):
        if position < (len(sequence) - 1):
            dinucleotide = character + sequence[position + 1]
            runningEnthalpy += meltingTable[dinucleotide]['H']
            runningEntropy += meltingTable[dinucleotide]['S']
    meltemp = (-(runningEnthalpy)/(-0.0108 - (runningEntropy/1000) + (0.00199 * math.log(oligoConcentration/4)))) - 273.15 + (16.6 * math.log10(sodiumConcentration))
    return meltemp

#as 3' and 5' junctions are sequenced differently, do not worry overmuch about disambiguating them. So focus on within each. Return frame offending sequences and number of times they appear
def uniquelyMapped(primerFrame, numberBases, ignoreBases):
    returnFrame = primerFrame.copy()
    returnFrame['basesRead'] = returnFrame.annealingSequence.str.slice(start=ignoreBases, stop=numberBases)
    primerFrame['uniqueInstances'] = primerFrame.groupby('annealingSequence').count().reset_index().segment
    primerFrame = primerFrame[primerFrame.uniqueInstances > 1]
    return primerFrame

#creates a concensus sequence between two reads using maximally scored positions
#returns concensus sequence, score, and the number of differences
def concensusSeq(seq1, seq2, score1, score2):
    concensus = ''
    concensusScore = ''
    difference = 0
    if len(seq1) == len(seq2):
        for position, character1, in enumerate(seq1):
            character2 = seq2[position]
            score1pos = score1[position]
            score2pos = score2[position]
            if character1 != character2:
                if score1pos > score2pos:
                    concensus += character1
                else:
                    concensus += character2
                difference += 1
            else:
                concensus += character1
            concensusScore += max((score1pos, score2pos))
    return (concensus, concensusScore, difference)







#hamming distance, ignoring N's in comparison
def seqHamm(sequence, comparison):
    hamm = 0
    for position, character in enumerate(sequence):
        if comparison[position] != 'N':
            if character != comparison[position]:
                hamm +=1
    return hamm
def seqHamQuick(sequence, comparison, limit):
    hamm = 0
    for position, character in enumerate(sequence):
        if comparison[position] != 'N':
            if character != comparison[position]:
                hamm +=1
        if hamm > limit:
            break
    return hamm

def minScore(scoreSeq):
    minVal = ord(scoreSeq[0]) - 33
    for character in scoreSeq:
        score = ord(character) - 33
        if score < minVal:
            minVal = score
    return minVal

def junctionAssign(*, Read1, Read2, UpMatch, DownMatch, outfile, Qvalcutoff=20, maxDistAdapter = 3,readDisagree = 1, junctionLen, barcodeLen):
    output = []
    with open(Read1, 'r') as read1File, open(Read2, 'r') as read2File:
        name1 = read1File.readline().split(' ')[0]
        name2 = read2File.readline().split(' ')[0]
        while len(name1) != 0:
            read1 = read1File.readline()[:-1]
            read2 = reverseComplement(read2File.readline()[:-1])
            toss = read1File.readline()
            toss = read2File.readline()
            score1 = read1File.readline()[:-1]
            score2 = read2File.readline()[:-1][::-1]
            upPos, upDist = minHam(UpMatch, read1)
            downPos, downDist = minHam(DownMatch, read1)
            #determine if an up or down flank, should be straightforward. Reannotate to make downstream code simpler

            if upDist < downDist:
                Up = True
                middlePos1 = upPos
                middlePos2, upDist2 = minHam(UpMatch, read2)
                middleDist = min(upDist, upDist2)
            else:
                Up = False
                middlePos1 = downPos
                middlePos2, downDist2 = minHam(DownMatch, read2)
                middleDist = min(downDist, downDist2)


            #junctions provided in primer orientation, opposite to adapter orientation
            junction1 = reverseComplement(read1[middlePos1-junctionLen:middlePos1])
            junction2 = reverseComplement(read2[middlePos2-junctionLen:middlePos2])
            scoreJunction1 = score1[middlePos1-junctionLen:middlePos1][::-1]
            scoreJunction2 = score2[middlePos2-junctionLen:middlePos2][::-1]
            junction, scoreJunction, junctionDisagree = concensusSeq(junction1,junction2, scoreJunction1,scoreJunction2)
            #if hamming distance for adapter appropriate
            if middleDist < maxDistAdapter:
                if Up:
                    #pull the barcode downstream of adapter
                    barcode1 = read1[middlePos1 + len(UpMatch):middlePos1 + len(UpMatch) + barcodeLen]
                    scoreBarcode1 = score1[middlePos1 + len(UpMatch):middlePos1 + len(UpMatch) + barcodeLen]
                    barcode2 = read2[middlePos2 + len(UpMatch):middlePos2 + len(UpMatch) + barcodeLen]
                    scoreBarcode2 = score2[middlePos2 + len(UpMatch):middlePos2 + len(UpMatch) + barcodeLen]
                    frame = pd.DataFrame({'orientation':['UP']})
                else:
                    #for Down flank, barcode is sequenced in opposite orientation, adjust
                    barcode1 = reverseComplement(read1[middlePos1 + len(DownMatch):middlePos1 + len(DownMatch) + barcodeLen])
                    scoreBarcode1 = score1[middlePos1 + len(DownMatch):middlePos1 + len(DownMatch) + barcodeLen][::-1]
                    barcode2 = reverseComplement(read2[middlePos2 + len(DownMatch):middlePos2 + len(DownMatch) + barcodeLen])
                    scoreBarcode2 = score2[middlePos2 + len(DownMatch):middlePos2 + len(DownMatch) + barcodeLen][::-1]
                    frame = pd.DataFrame({'orientation':['DOWN']})
                barcode, scoreBarcode, barcodeDisagree = concensusSeq(barcode1,barcode2, scoreBarcode1,scoreBarcode2)
                if (len(barcode) == barcodeLen) and (len(junction) == junctionLen):
                    barcodeScoreMinimum = minScore(scoreBarcode)
                    junctionScoreMinimum = minScore(scoreJunction)
                    if (junctionScoreMinimum >= Qvalcutoff) and (barcodeScoreMinimum >= Qvalcutoff) and (junctionDisagree <= readDisagree) and (barcodeDisagree <= readDisagree):
                        frame['junction'] = junction
                        frame['barcode'] = barcode
                        frame['count'] = 1
                        output += [frame]
            name1 = read1File.readline().split(' ')[0]
            name2 = read2File.readline().split(' ')[0]

    pd.concat(output).groupby(['junction', 'orientation', 'barcode']).sum().to_csv(outfile, sep='\t')

#collapse junctions to set of known junctions. Allow a certain hamming distance. Require matches UP or DOWN as indicated by sequencing
def collapseJunctions(*, junctions, infile, permittedDistance, outfile):
    outputFrame = []
    UPjunctions = junctions.loc[junctions.orientation == 'UP',]
    DOWNjunctions = junctions.loc[junctions.orientation == 'DOWN']
    barcodeFrame = pd.read_csv(infile, sep='\t')
    for index, row in barcodeFrame.iterrows():
        junction = row.junction
        distance = len(junction)
        if row.orientation == 'UP':
            relevant = UPjunctions
        else:
            relevant = DOWNjunctions
        if len(relevant.loc[relevant.sequenceInRead == junction,]) != 1:
            for indexRow, junctionRow in relevant.iterrows():
                currJunction = junctionRow.sequenceInRead
                currDist = seqHamQuick(junction,currJunction, permittedDistance)
                if currDist < distance:
                    distance = currDist
                    correctedJunction = indexRow
            if distance <= permittedDistance:
                outFrame = pd.DataFrame(relevant.loc[correctedJunction])
                outFrame['count'] = row['count']
                outFrame['barcode'] = row.barcode
                outputFrame += [outFrame]

        else:
            outFrame = pd.DataFrame(relevant.loc[relevant.sequenceInRead == junction,])
            outFrame['count'] = row['count']
            outFrame['barcode'] = row.barcode
            outputFrame += [outFrame]
    outputFrame = pd.concat(outputFrame)
    #collapse identical barcodes that had errors in junctions
    outputFrame = outputFrame.groupby(['barcode', 'segment','orientation', 'position', 'sequenceInRead']).sum().reset_index()
    outputFrame.to_csv(outfile, sep='\t')




#output raw barcode counts, can do processing later. Barcode is in forward orientation for R2. Require both reads match barcode. Use minimum score on either read to define threshold.
#changed from original to not use pd concat as I think it was giving a memory issue, instead use a library collapsed to pandas at end.
def barcodeCount(*, Read1, Read2, adjSeq, Qvalcutoff=20, maxDistAdapter = 3, barcodeLen,  readDisagree = 1, outfile):
    data = {}
    with open(Read1) as R1, open(Read2) as R2:
        line1 = R1.readline()
        line2 = R2.readline()
        while line1 != '' and line2 != '':
            #files need to have equal numbers of reads, so check that out
            try:
                assert line1.split(' ')[0] == line2.split(' ')[0]
                line1 = R1.readline()[:-1]
                line2 = R2.readline()[:-1]
                line1 = reverseComplement(line1)
                score1 = R1.readline()
                score1 = R1.readline()[:-1][::-1]
                score2 = R2.readline()
                score2 = R2.readline()[:-1]
                barcodePos1, mismatch1 = minHam(adjSeq, line1)
                barcodePos2, mismatch2 = minHam(adjSeq, line2)
                if mismatch1 <= maxDistAdapter and mismatch2 <= maxDistAdapter:
                    barcode1 = line1[barcodePos1 + len(adjSeq): barcodePos1 + len(adjSeq) + barcodeLen]
                    score1 = score1[barcodePos1 + len(adjSeq): barcodePos1 + len(adjSeq) + barcodeLen]
                    barcode2 = line2[barcodePos2 + len(adjSeq): barcodePos2 + len(adjSeq) + barcodeLen]
                    score2 = score2[barcodePos2 + len(adjSeq): barcodePos2 + len(adjSeq) + barcodeLen]
                    #pull concensus, mismatch, score
                    barcode, scoreBarcode, barcodeDisagree = concensusSeq(barcode1,barcode2, score1,score2)
                    if barcodeDisagree <= readDisagree:
                        if len(barcode) == barcodeLen:
                            if minScore(scoreBarcode) >= Qvalcutoff:
                                currData = data
                                for letter in barcode:
                                    if letter not in currData:
                                        currData[letter] = {}
                                    currData = currData[letter]
                                if len(currData) == 0:
                                    currData['data'] = pd.DataFrame({'barcode':[barcode], 'count':[1]})
                                else:
                                    currData['data']['count'] += 1

            except AssertionError:
                print('Read names discordant. Offending names are ' + line1 + ' and ' + line2 + '. R1 and R2 files must be sorted by read name and all be paired')
            line1 = R1.readline()
            line2 = R2.readline()
    output = []

    data = pd.concat(NestedDictValues(data))
    data.to_csv(outfile, sep='\t', index=False)




#returns the 0 position of the best match and the hamming distance as a tuple
def minHam(target, search):
    currPos = 0
    minPos = len(target)
    minHam = len(target)
    while(len(search) - currPos >=  len(target)):
        ham = seqHamm(target, search[currPos:currPos+len(target)])
        if ham < minHam:
            minHam = ham
            minPos = currPos
        currPos += 1
    return (minPos, minHam)

#from https://stackoverflow.com/questions/23981553/get-all-values-from-nested-dictionaries-in-python
def NestedDictValues(d):
  for v in d.values():
    if isinstance(v, dict):
      yield from NestedDictValues(v)
    else:
      yield v


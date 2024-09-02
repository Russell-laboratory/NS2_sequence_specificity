
import pandas as pd
import os

#list of bases recorded at each position that meet a minimum qscore. can then compare against genome later. 1-index rather than 0
def positionalData(*, sortedBamfile, outfile, qvalCutoff):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 2 -F 524', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                insertions ={}
                segment = firstMate[2]
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            #only count if quality score to either side meets standard. Only called if base not already called at position.
                            if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                                while basesAdded < value[0]:
                                    currentSeq[currentPosSeq] = ('D', qvalCutoff)
                                    basesAdded += 1
                                    currentPosSeq += 1
                                    #do not advance place in read
                            else:
                                currentPosSeq += value[0]
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        #inserted sequence gets added at the position if bases before and after are a good score
                        elif value[1] == 'I':
                            if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                                insertions[currentPosSeq] = value[0]
                                currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]

                for position in currentSeq:
                    if segment not in output:
                        output[segment] = {}
                    if position not in output[segment]:
                        output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                    if currentSeq[position][1] >= qvalCutoff:
                        #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                        base = currentSeq[position][0]
                        if base != 'N':
                            output[segment][position][base] += 1
                for position in insertions:
                    output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"

#ignore positions covered by primers. Use known subamplicon coordinates and ignore given positions. Ignorepos given as a dictionary.
def positionalDataTrim(*, sortedBamfile, outfile, qvalCutoff, ignorePos):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 2 -F 524', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                insertions ={}
                segment = firstMate[2]
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            #only count if quality score to either side meets standard. Only called if base not already called at position.
                            if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                                while basesAdded < value[0]:
                                    currentSeq[currentPosSeq] = ('D', qvalCutoff)
                                    basesAdded += 1
                                    currentPosSeq += 1
                                    #do not advance place in read
                            else:
                                currentPosSeq += value[0]
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        #inserted sequence gets added at the position if bases before and after are a good score
                        elif value[1] == 'I':
                            if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                                insertions[currentPosSeq] = value[0]
                                currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]
                ignoreList = []
                #now "trim" Use a tuple-keyed list of positions unique to each amplicon and then ignore positions defined in thing
                for amplicon in ignorePos.keys():
                    if currentSeq.keys() & amplicon:
                        ignoreList += ignorePos[amplicon]
                if len(ignoreList) != 0:
                    for position in currentSeq:
                        if position not in ignoreList:
                            if segment not in output:
                                output[segment] = {}
                            if position not in output[segment]:
                                output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                            if currentSeq[position][1] >= qvalCutoff:
                                #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                                base = currentSeq[position][0]
                                if base != 'N':
                                    output[segment][position][base] += 1
                    for position in insertions:
                        if position not in ignoreList:
                            output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"

#quick and dirty to handle caps before I write something better that takes it pre-alignemnt
def positionalDataIgnoreInsert(*, sortedBamfile, outfile, qvalCutoff, insertionCutoff):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 2 -F 524', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                insertions ={}
                segment = firstMate[2]
                insertionTooLarge = False
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            if character == 'I':
                                if int(curstring) >= insertionCutoff:
                                    insertionTooLarge = True
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            #only count if quality score to either side meets standard. Only called if base not already called at position.
                            if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                                while basesAdded < value[0]:
                                    currentSeq[currentPosSeq] = ('D', qvalCutoff)
                                    basesAdded += 1
                                    currentPosSeq += 1
                                    #do not advance place in read
                            else:
                                currentPosSeq += value[0]
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        #inserted sequence gets added at the position if bases before and after are a good score
                        elif value[1] == 'I':
                            if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                                insertions[currentPosSeq] = value[0]
                                currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]
                ignoreList = []
                #now "trim" Use a tuple-keyed list of positions unique to each amplicon and then ignore positions defined in thing
                if not insertionTooLarge:
                    for position in currentSeq:
                        if position not in ignoreList:
                            if segment not in output:
                                output[segment] = {}
                            if position not in output[segment]:
                                output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                            if currentSeq[position][1] >= qvalCutoff:
                                #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                                base = currentSeq[position][0]
                                if base != 'N':
                                    output[segment][position][base] += 1
                    for position in insertions:
                        if position not in ignoreList:
                            output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"

#assumes sequence of interest is in R2, if not, invert
def cappedExclude(*, searchSeqs, R1in, R2in, R1out,R2out, toleratedDifs):
    with open(R1in) as R1in, open(R2in) as R2in, open(R1out, 'w') as R1out, open(R2out, 'w') as R2out:
        name1 = R1in.readline()
        name2 = R2in.readline()
        while name1 != '':
            seq1 = R1in.readline()
            seq2 = R2in.readline()
            plus1 = R1in.readline()
            plus2 = R2in.readline()
            score1 = R1in.readline()
            score2 = R2in.readline()
            #find best match. Does not work for indels, oh well. Repeat for plasmid samples for parity then.
            minDifs = toleratedDifs + 1
            for sequence in searchSeqs:
                position = len(sequence)
                sub = seq2[position - len(sequence):position]
                #borrowed from https://stackoverflow.com/questions/28423448/counting-differences-between-two-strings
                difs = sum(1 for a, b in zip(sequence, sub) if a != b)
                if minDifs > difs:
                    minDifs = difs
            if minDifs <= toleratedDifs:
                R1out.write(name1)
                R1out.write(seq1)
                R1out.write(plus1)
                R1out.write(score1)
                R2out.write(name2)
                R2out.write(seq2)
                R2out.write(plus2)
                R2out.write(score2)
            name1 = R1in.readline()
            name2 = R2in.readline()
                    

#list of bases recorded at each position. Assume from PrimerID that you can safely ignore qVal. can then compare against genome later. 1-index rather than 0
def positionalDataSE(*, bamfile, outfile):
    tempFile = bamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -F 8', bamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        for line in infile:
            read = line.split('\t')            
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            currentSeq = {}
            insertions ={}
            segment = read[2]
            CIGAR = read[5]
            currentPosSeq = int(read[3])
            score = read[10]
            sequence = read[9]
            #start positions and mapping.
            CIGARval = []
            curstring = ''
            for character in CIGAR:
                if not character.isalpha():
                    curstring += character
                else:
                    CIGARval.append((int(curstring), character))
                    curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                currentSeq[currentPosSeq] = sequence[currentPosRead-1]
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            while basesAdded < value[0]:
                                currentSeq[currentPosSeq] = 'D'
                                basesAdded += 1
                                currentPosSeq += 1
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        elif value[1] == 'I':
                            insertions[currentPosSeq] = value[0]
                            currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]

                for position in currentSeq:
                    if segment not in output:
                        output[segment] = {}
                    if position not in output[segment]:
                        output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                    base = currentSeq[position]
                    if base != 'N':
                        output[segment][position][base] += 1
                for position in insertions:
                    if segment not in output:
                        output[segment] = {}
                    if position not in output[segment]:
                        output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                    output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"


def fastaReport(*, inFasta,outfile):
    #support multiline fasta, although I wish that wasn't a thing
    with open(inFasta) as infile, open(outfile, 'w') as outfile:
        positions = {}
        currSeq = ''
        for line in infile:
            line = line[:-1]
            if line[0] == '>':
                for position, character in enumerate(currSeq):
                    position = position + 1
                    if position not in positions:
                        positions[position] = {'A':0, 'G':0,'C':0,'T':0}
                    positions[position][character] += 1
                currSeq = ''
            else:
                currSeq += line
        for position, character in enumerate(currSeq):
            position = position + 1
            if position not in positions:
                positions[position] = {'A':0, 'G':0,'C':0,'T':0}
            positions[position][character] += 1
        outfile.write('\t'.join(['position','A','G','C','T']) + '\n')
        for position in positions:
            curr = positions[position]
            outfile.write('\t'.join([str(position), str(curr['A']), str(curr['G']), str(curr['C']), str(curr['T'])]) + '\n')
    



        
#scripts rely heavily on fastQC being installed in samtools_path

import subprocess
import pandas as pd
from io import StringIO
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

'''runs fastQC on both given set of sample-associated reads. Returns
datasets for graphing. Reads provided as a dictionary with names'''
def fastQCreport(*,sampleName, reads):
    #check if file exists, if it does we can just grab the data instead
    for name in reads.keys():
        read = reads[name]
        if not os.path.isfile(read.split('.')[0] + '_fastqc/fastqc_data.txt'):
            process = subprocess.Popen(' '.join(['fastqc ', read, '--extract']), shell=True)
            process.communicate()
    #now just grab the stuff you need from the files
    pandasReports = []
    x = 1
    for name in reads.keys():
        read = reads[name]
        report = {}
        statistics={}
        with open(read.split('.')[0] + '_fastqc/fastqc_data.txt', 'r') as infile:
            #get rid of hashtags as they are annoying as all get out
            content = infile.read().replace("#", "").split('>>')
        report['file'] = str(x)
        x +=1
        for item in content:
            try:
                #title of section, remove pass statistic
                header = item.split('\n')[0].split('\t')[0]
                #meat of section
                data = pd.read_csv(StringIO('\n'.join(item.split('\n')[1:])), sep='\t')
                statistics[header] = data
            except Exception:
                pass
        #number of sequences alone is fine
        basicStatistics = statistics['Basic Statistics']
        report['totalSeq'] = int(basicStatistics[basicStatistics.Measure == 'Total Sequences'].Value.iloc[0])
        #now length distribution
        report['lengthDist'] = statistics['Sequence Length Distribution']
        # a bit of sillyness due to mixed types, cast it as a string only to cast it back later
        report['lengthDist']['Length'] = report['lengthDist']['Length'].astype("str")
        #when report is split across two bases, just take first one. It honestly won't make a difference and makes this mappable.
        report['lengthDist']['Length']=report['lengthDist'].apply(lambda row: int(row['Length']) if row['Length'].isdigit() else int(row['Length'].split('-')[0]), axis=1)

        #QC Distribution
        report['QClength'] = statistics['Per base sequence quality'].melt(id_vars=['Base'], value_vars=['Mean']).rename(columns={'value':'mean Qscore'})
        #Same sillyness as above
        report['QClength']['Base'] = report['QClength'].Base.astype("str")
        report['QClength']['Base'] = report['QClength'].apply(lambda row: int(row.Base) if row.Base.isdigit() else int(row.Base.split('-')[0]), axis=1)
        #now QC per read distribution
        report['QCread'] =    statistics['Per sequence quality scores']
        report['readName'] = name
        pandasReports += [report]
    return pandasReports

'''Handler to generate a graphical summary of the state of the reads used in this study from fastQC output.
A summary barplot of all samples, all reads, of just amount of sequencing. Then break
out for all other stuff. each sample alloted a row, then subdivided by readsets.'''
def graphAllTheThings(*,data, sampleFont=16, color='0.5',
depthHeight=2, columnwidth=4, lengthdistHeight=2, lengthdistsize=50, lengthQCheight=3, lengthQCsize=50, QCreadHeight=3, QCreadSize=50,
titleHeight = 0, titleWidth = 4, title='', titleFont = 28):
    #rows comprising a row. Set to max reads in a sample. WIll just leave blank any additional areas
    rowSub = max([len(data[key]) for key in data.keys()])
    #3 categories with subdivided rows, depth and a title row
    rows = rowSub*3 + 2
    columns = len(data.keys()) + 1
    #initialize figure, then gridspec using the appropriate heights and widths
    g= plt.figure(1,((columns * columnwidth), (lengthdistHeight  + lengthQCheight  + QCreadHeight) * rowSub + depthHeight + titleHeight))
    widths = [titleWidth] + [columnwidth] * (columns - 1)
    heights = [titleHeight, depthHeight] + [lengthdistHeight] * rowSub + [lengthQCheight] * rowSub + [QCreadHeight] * rowSub
    spec = g.add_gridspec(ncols=columns,nrows=rows, width_ratios = widths, height_ratios=heights)

    #lets get those titles in there. First columns
    columnTitles = data.keys()
    col = 1
    for sample in columnTitles:
        ax =  g.add_subplot(spec[0, col])
        ax.axis('off')
        ax.text(0,-0.5, sample, transform=ax.transAxes, fontsize=sampleFont)
        col += 1
    #now lets get the row titles in there. For subdivided rows, summate rows to center title
    ax = g.add_subplot(spec[1, 0])
    ax.axis('off')
    ax.text(0,0, 'read count', transform=ax.transAxes, fontsize=sampleFont)
    subdivRowOrder = ['length distribution', 'QC over read length', 'QC distribution']
    row = 2
    for element in subdivRowOrder:
        ax = g.add_subplot(spec[row:row+rowSub, 0])
        ax.axis('off')
        ax.text(0,0.5, element, transform=ax.transAxes, fontsize=sampleFont)
        row += rowSub

    row=1
    col=1
    topAxes = {}
    for sample in data.keys():
        #first just bargraph of number of reads for all samples in samples
        readLengths = pd.DataFrame(columns=['file','reads'])
        for read in data[sample]:
            readLengths = pd.concat([readLengths, pd.DataFrame(data = {'file':[read['readName']], 'reads':[read['totalSeq']]})])
        if 'numSeq' not in topAxes.keys():
            ax = g.add_subplot(spec[row, col])
            topAxes['numSeq'] = ax
        else:
            ax = g.add_subplot(spec[row, col], sharey=topAxes['numSeq'])
        sns.barplot(data=readLengths, x='file', y='reads', color=color, ax=ax)
        ax.tick_params(axis='x', rotation=45)
        ax.yaxis.label.set_size(sampleFont)
        ax.xaxis.label.set_size(sampleFont)
        col += 1

    #very repetitive as rest are all scatterplots by read, so go ahead and store everything that varies here and just adjust
    stuffToIterate = (('lengthDist',('Length', 'Count', lengthdistsize)),
                        ('QClength',('Base','mean Qscore',lengthdistsize)),
                        ('QCread', ('Quality','Count', QCreadSize )))
    col = 1
    for sample in data.keys():
        row = 2
        for metric in stuffToIterate:
            numRowsAdded = 0
            for read in data[sample]:
                if metric[0] not in topAxes.keys():
                    ax = g.add_subplot(spec[row + numRowsAdded,col])
                    topAxes[metric[0]] = ax
                else:
                    ax = g.add_subplot(spec[row + numRowsAdded, col], sharex=topAxes[metric[0]], sharey=topAxes[metric[0]])
                numRowsAdded += 1
                sns.scatterplot(data=read[metric[0]], x=metric[1][0], y=metric[1][1], color=color, ax=ax,  s=metric[1][2])
                ax.yaxis.label.set_size(sampleFont)
                ax.xaxis.label.set_size(sampleFont)
            #increment row the maximum number of reads regardless of how many were in this sample
            row += rowSub
        #next sample
        col += 1
    plt.tight_layout()
    if len(title) != 0:
        g.suptitle(title, fontsize=titleFont, y=1.1)
    return g

#simple handler to call trimmomatic. If files already exist does not overwrite. Returns statistics from trimmomatic. Writes statistics to file for subsequent reruns.
def trimCommandPaired(*, read1, read2, outDirectory, adapterFile, seedMismatches=2, palindromeClipThreshold=30, simpleClipThreshold=10,
minAdapterLength=2, keepBoth="TRUE", lead=20, slidingMin=4, slidingMax=15, threads=1, minlentrim=36, overwrite=False, sampleName = ''):

    #use pretty standard R1 and R2 designations. If becomes a problem later can always pass as a variable
    read1pairedtrim = outDirectory + '/' + read1.split('/')[-1]
    read1unpairedtrim= outDirectory + '/Unpaired_' + read1.split('/')[-1]
    read2pariedtrim = outDirectory + '/' + read2.split('/')[-1]
    read2unpairedtrim = outDirectory + '/Unpaired_' + read2.split('/')[-1]
    terminalCapture = outDirectory + '_' + read1.split('/')[-1].replace('_R1','').replace('_R2','').split('.')[0] +  "_logfile.txt"
    #just check for paired trimmed files, hopefully not to much in the unpaired files if I am doing this correctly. If there is..definitely need to recheck parameters
    if not ((os.path.isfile(read1pairedtrim) & os.path.isfile(read2pariedtrim)) or overwrite):
        illuminaClipArgument = ':'.join(['ILLUMINACLIP',adapterFile, str(seedMismatches), str(palindromeClipThreshold), str(simpleClipThreshold), str(minAdapterLength),keepBoth])
        arguments = ' '.join(["trimmomatic PE -threads", str(threads),
            '-phred33', read1, read2, read1pairedtrim, read1unpairedtrim,
                read2pariedtrim, read2unpairedtrim, illuminaClipArgument, "LEADING:" + str(lead),
                "SLIDINGWINDOW:" + str(slidingMin) + ':' + str(slidingMax), "MINLEN:" + str(minlentrim)])
        #a bit bad of me but I am just using shell=True here as subprocess for java jar was giving me a headache. This isn't general software so a bit less
        #concerned about shell injection!
        process = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stderr, stdout = process.communicate()
        stdout = stdout.decode('utf-8')
        #capture any issues for posterity
        with open(terminalCapture, 'w') as out:
            out.write(stdout)

    with open(terminalCapture, 'r') as infile:
        stats = infile.readlines()[-2]
        stats = stats.split(' ')
    return pd.DataFrame(data={'Sample':sampleName, 'Input read pairs':[int(stats[3])], 'Both surviving':[int(stats[6])], 'Forward only':[int(stats[11])],
            'Reverse only':[int(stats[16])], 'dropped':[int(stats[-2])]})


#returns dataframe
def STARreport(infile, sampleName):

    data = pd.DataFrame()
    data['Sample name'] = [sampleName]
    with open(infile, 'r') as infile:
        for line in infile:
            line = line.strip()
            lineseg = line.split('|')
            if len(lineseg) == 2:
                if lineseg[1][-1] != '%':
                    data[lineseg[0].strip()] = [lineseg[1].strip()]
                else:
                    data[lineseg[0].strip()] = [float(lineseg[1].strip()[:-1])]

    return data

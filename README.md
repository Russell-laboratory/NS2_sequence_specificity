

# NS2_sequence_specificity

The code contained within this repository describe the analysis of data, and generation of figures, in our attempts to understand how expression of NS2 influences the expression of libraries of length-variant HA and PB1, and sequence-variant PB1.

Analysis are split across several different jupyter notebooks, with descriptions below, and figures inline.

All intermediate files less than 10mb were pushed to this repository to aid in understanding our analysis.
Larger files, such as sequencing, will need to be downloaded and reconstructed in the apppropriate folders in order to re-run this analysis.
Specifically, place NGS files within appropriate directories within Sequencing to rerun all sequencing analyses.

Data may be found in GEO (including some processed files) at accession GSE276697.

## Directories

The following directories, and their purpose, exist within this repository.

- <b>Database</b>       Repository for influenza genomic sequences and the same with adapters for appropriate mapping as well as sequence with plasmid backbone. A/WSN/1933 BLAST database and STAR indicies provided.
-  <b>Figures</b>       Figures, main and some supplemental.
- <b>SuppFigures</b>    Additional supplemental figures
- <b>Results</b>        Final datafiles for analyses after processing. Most critical ones additionally on GEO.
- <b>Scripts</b>        Short scripts written for this analysis. Seperated from jupyter notebooks for readability and portability.
- <b>Sequencing</b>     Folder that would contain NGS samples. Folder achitecture essential to run off pipeline.
- <b>PriorData</b>      Prior data from Mendes et al. 2021 for analyzing length polymorphic libraries (identical libraries were used as generated in this paper).
  

## Dependencies

Code within this repository was run with the following tools, and versions, installed and available from PATH. Specific websites and documentation are provided where available. 

- <b>Python</b>      run with version 3.12.2 . Available from https://www.python.org/downloads/release/python-370/
- <b>Trimmomatic</b> run with version 0.39. Available from http://www.usadellab.org/cms/?page=trimmomatic
- <b>STAR</b>        run with version 2.7.11b. Available from https://github.com/alexdobin/STAR
- <b>Samtools</b>    run with version 1.19. Available from http://www.htslib.org/
- <b>FastQC</b>      run with version 0.12.1. Available from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


The following python packages and versions were used. All were installed using conda. (https://docs.conda.io/en/latest/)
- <b>numpy</b>       run with version 1.26.4. (https://numpy.org/)
- <b>matplotlib</b>  run with version 3.8.3. (https://matplotlib.org/)
- <b>seaborn</b>     run with version 0.13.2. (https://seaborn.pydata.org/)
- <b>pandas</b>      run with version 2.2.1.(https://pandas.pydata.org/)
- <b>scipy</b>       run with version 1.12.0. (https://www.scipy.org/)
- <b>statsmodels</b> run with version 0.14.1. (https://www.statsmodels.org/stable/index.html)

## Jupyter notebooks

General descriptions of pipelines within each notebook described below.

### Figures_sequencing.ipynb 

Generating figures and interpreting raw data from other notebooks.

### Length_polymorphic_assignation.ipynb

Extraction of barcodes from sequencing data of length-variant libraries. Uses code identical to Mendes et al. 2021.

### Natural_diverisity
Analysis of diverse PB1 sequences from the NCBI flu database. Preprocessing performed outside of this notebook included CD-HIT and manual curation.


### SNP mapping and enumeration
Alignment to either PB1_177_385, the same with a TSO (for RACE), or the same with plasmid backbone (RACE plasmid comparison).
For non-RACE samples sequences processed with trimmomatic prior to alignment to appropriate template with STAR.
For non-RACE samples, used a quality score cutoff of 30 and enumerated all nucleotides observed at a given position.
For RACE-samples, reads were processed by removing those with greater than 4 (plasmid comparison) or 5 (RACE) mismatches to the anticipated 5' sequence of either vRNA or cRNA. 
This step was critical to remove cap-containing sequences. 
After this curation, STAR was used to map remaining sequences, and nucleotides at each position called with a q-score cutoff of 30 (with additional exclusion of any remaining mapped reads that contained large, heterogeneous, sequences consistent with cap). 



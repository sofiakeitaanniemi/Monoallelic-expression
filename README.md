# Monoallelic-expression
Bachelor's thesis project: Python analysis workflow for monoallelic expression analysis

## In short
Samtools pileup file parsing, combining parsed data with allele data and finding SNPs located in genes. 

## HISTORY
In an example analysis, Affymetrix SNP 6.0 array data (CEL files) were downloaded and 
processed with R package Bioconductor Crlmm (code included). As output we got allele calls for SNPs. 
We also needed coverage data from GRO-SEQ and RNA-SEQ BAM files. BAM files were run with Samtools mpileup and 
the output PILEUP files were run with CountReads_fromPileup.py. Allele calls and 
parsed pileup data were combined in AlleleFinding_forParsedPileup.py. 
Finally, we looked for SNPs in genes: we counted the number of SNPs in genes (CSV output file) 
and their positions (BED) using SNP_inGenes.py. Reference gene file was downloaded from 
UCSC Gene Table as BED file (RefSeq Genes).

Created with Python 3.6
13/07/2017 


## Sub-modules
CountReads_fromPileup.py, AlleleFinding_forParsedPileup.py, SNP_inGenes.py
Also included: Bioconductor_crlmm R file for retrieving SNP alleles (calls) from Affymetrix SNP 6.0 array CEL files.
The code creates tab delimited text file. 
Sub-modules are supposed to be used in the following order. 
If only some parts of the project are used, take into account the input explanations 
because some of the needed files are created in the pervious module. 
Check that the files you are using in your analysis are suitable for input.


### 1. CountReads_fromPileup.py
The first module reads an input Samtools pileup file created from GRO-SEQ or RNA-SEQ BAM files 
(usage: samtools mpileup -l SNPfile.bed -f reference_genome.fa Gro/Rna-seq_BAMfile.bam > Your_new_PILEUPfile.pileup).
The module counts the reads and creates a new tab delimited file (BED). The file name is asked from the user and
extension '_countedPileup.bed' is added to the end of the file name. Caution: The output file cannot be
viewed in IGV. For that you need to run the next python file or change the code and in that way the number of columns
written to this first BED file.
Needed packages: Counter from collection and re which are imported in python file.

USAGE:
python3 CountReads_fromPileup.py

INPUT EXPLANATION: 
Samtools pileup file

OUTPUT EXPLANATION:
Tab delimited file (BED)
"CHR START END REF NUM A T C G N %-A %-T %-C %-G %-N"



### 2. AlleleFinding_forParsedPileup.py
The second module needs a tab delimited file as input (created in CountReads_fromPileup.py or similar). 
The second input, SNP Alleles, is created with R package Bioconductor Crlmm (example code in this 
MonoallelicExpression package) using Affymetrix SNP 6.0 array CEL files, for example. 
The input file must be tab delimited (BED), consisting of data: 
chromosome, SNP start position, SNP end position and allele (1/2/3).
Output file is tab delimited BED file. The file name is asked from the user and 
an extension '_alleleCoverage.bed' is added to the end of the file name. 
The number of columns in newly created file is lower than in the first module output
so that the files can be viewed in IGV.

USAGE:
python AlleleFinding_forParsedPileup.py

INPUT EXPLANATION:
 1. Parsed pileup file: tab delimited BED file created with CountReads_fromPileup.py
"CHR START END REF NUM A T C G N %-A %-T %-C %-G %-N"
 2. SNP Alleles: tab delimited BED file (created from CEL files with R package Bioconductor Crlmm, for example)
"CHR START END ALLELE"
3. Select which alleles are processed (1, 2, 3 or all). 1 = AA. 2 = AB, 3 = BB
4. Enter percentage treshold after allele selection. Only SNPs where a singular nucleotide has percentage 
of reads over treshold, are saved for further analysis. Select zero (0) if you want to analyze all the SNPs.  
5. Enter name for the new BED file

OUTPUT EXPLANATION:
Tab delimited file (BED)
"CHR START END REF NUM %-A %-T %-C %-G %-N ALLELE"



### 3. SNP_inGenes.py
The final module reads gene reference file and the tab delimited BED file created previous python file, 
AlleleFinding_forParsedPileup.py. Two files are created. One is a CSV file of gene names and number of SNPs
located in these genes. The other is a tab delimited BED file consisting of chromosome, start of the SNP, 
end of the SNP and gene name.
Gene reference file must be tab delimited (BED) file with GENE ID, CHR, START, END and GENE NAME. 
This kind of file can be downloaded from UCSC Gene Table, for example. 

USAGE:
python SNP_inGenes.py

INPUT EXPLANATION:
1. Allele coverage file: tab delimited BED file created with AlleleFinding_forParsedPileup.py
"CHR START END REF NUM %-A %-T %-C %-G %-N ALLELE"
2. Gene reference file: tab delimited BED file (USCS Gene Table, RefSeq genes, for example)
"REFSEQ_ID CHR START END GENE_NAME"
3. Enter name for the new CSV file. CSV file consists of gene names and the number of SNPs in each gene.
4. Enter name for the new BED position file. 

OUTPUT EXPLANATION:
1.(file name from 3. input) CSV file consist of gene names and the number of SNPs in each gene.
"GENE_NAME,NUMBER_OF_SNPs"
2.(file name from 4. input) BED file consist of Gene names and positions for each SNP
"CHR START END Gene_name"

In the end of the last run there are two files created: SNP positions and number of SNPs in genes
 depending on the selected allele and percentage treshold. 



## Author / Contact:
Sofia Randelin, University of Tampere, Computational biology
sofia.randelin@tuni.fi

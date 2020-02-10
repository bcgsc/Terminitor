# Terminitor

Terminitor is a deep neural network that predicts whether a sequence contains a polyadenylated (poly(A)) cleavage site (CS) at certain position.

For more information, please refer to the preprint: https://www.biorxiv.org/content/10.1101/710699v1

### Datasets for download  
www.bcgsc.ca/downloads/supplementary/Terminitor  

This ftp site contains two datasets, human and mouse, and two corresponding pre-trained models for test.  

### Dependencies  
* Python3
* Numpy
* Keras
* Scikit-learn  
* Pybedtools  
* Pysam  

### Train  

Usage: `train.py [-h] [-v] -polya POLYA -cs CS -non NON -model MODEL -l L`

```
optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
  -polya POLYA   Poly(A) CS, fasta file
  -cs CS         Non-poly(A) CS, fasta file
  -non NON       Non-CS, fasta file
  -model MODEL   File name of trained model
  -l L           Length of input sequences
```

### Test

Usage: `test.py [-h] [-v] -t TEST_FILE -m MODEL -l L -o OUTPUT`  

```
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -t TEST_FILE, --test_file TEST_FILE
                        Fasta file to be tested
  -m MODEL, --model MODEL
                        Pre-trained model file
  -l L                  Length of input sequences
  -o OUTPUT, --output OUTPUT
                        Output probabilities
```

### Pipeline

1. For Illumina RNA-seq short reads, run assembly with [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom)  
`java -jar RNA-Bloom.jar -left read2.fq -right read1.fq -revcomp-right -outdir assembly -a 4 -e 1 -stratum 01 -ss -ntcard -fpr 0.005`  
For PacBio CCS reads, skip this step  

2. Genome alignment with [minimap2](https://github.com/lh3/minimap2)    
`minimap2 -ax splice hg38.mmi rnabloom.transcripts.fa | samtools view -b - > aln.bam`  

3. Extract candidate sequence  
`./src/extract_from_sequences.py -t Homo_sapiens.GRCh38.94.transcripts.gtf -a Homo_sapiens.GRCh38.94.chr.gtf -g Homo_sapiens.GRCh38.dna.primary_assembly.fa -m aln.bam -o test`  
Details are explained below.  

4. Test  
`./test.py -t extracted_sequence.fa -m pre_trained_model -l 200 -o output`  

### Extract candidate sequence  

Usage: `extract_from_sequences.py [-h] [-v] -t ANNOT_TRANS -a ANNOT_ALL -m ALN
                                 -g GENOME -o O [-u UP_LEN] [-d DOWN_LEN]`  

```
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -t ANNOT_TRANS, --annot_trans ANNOT_TRANS
                        Transcript annotation file, GTF format. This file
                        contains onlytranscript level annotation, can be
                        downloaded from the ftp siteprovided on our Github
                        page
  -a ANNOT_ALL, --annot_all ANNOT_ALL
                        Ensembl annotation file, GTF format. Can be downloaded
                        from Ensemblftp site
  -m ALN, --aln ALN     The alignment file from assembled transcript contigs
                        to reference genomein bam format.
  -g GENOME, --genome GENOME
                        Reference genome assembly in Fasta format. Can be
                        downloaded fromEnsembl ftp site
  -o O                  Output file, fasta format containing candidate
                        sequences to be tested
  -u UP_LEN, --up_len UP_LEN
                        Upstream sequence length
  -d DOWN_LEN, --down_len DOWN_LEN
                        Downstream sequence length
```

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


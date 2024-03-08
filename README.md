# MOTIF MARK

### INTRODUCTION

Motifs are biological sequences that are target for DNA/protein binding in regulatory elements. It can be helpful to visualize where different motifs are in a sequence. 

motif-mark-oop.py utilizes object oriented programming to take in a FASTA file (with wrapped DNA sequence lines), along with a text file containing all motifs of interest to generate an image to visualize exons and different motifs on the sequence. 

This program recognizes ambiguous nucleotide coding such as R,V,Y etc. in addition to A,T,G,C,U bases in motifs. List of other IUPAC degenerate base symbols can be found here: https://en.wikipedia.org/wiki/Nucleic_acid_notation


### GETTNG STARTED

INSTALLATION REQUIREMENTS:

Create a new environment: ```conda create -n my_pycairo pycairo```
 
Install pycairo: ```conda install pycairo```

Dependencies: Python version 3 or higher


HOW TO RUN:

```
conda activate my_pycairo

./motif-mark-oop.py -f <path to FASTA file> -m <path to motifs text file>
```

INPUT:

FASTA file (seqs ≤1000 bases) and 
motifs file (≤10 bases each, one motif per line in a text file)
- capable of handling maximum of 5 motifs




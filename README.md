# HGTmeta
An efficient algorithm for the inference of the optimal gene-species mappings under the under the DTL model.

## Requirements:
Python 2.7

## Usage:
python hgtmeta.py -S \<file with a species tree> -G \<file with a gene tree> [-D \<int>, -L \<int>, -T \<int>, -h]
  
###### Examples:
<br/>
For options description:

python hgtmeta.py -h

<br/>
Running the program with default options:
python hgtmeta.py -G data/Gtree_A -S data/Stree_A

<br/>
Setting costs for evolutionary events (duplication cost = 2, loss cost = 1, hgt cost = 3):
python hgtboot.py -G data/Gtree_A -S data/Stree_A -D 2 -L 1 -T 3

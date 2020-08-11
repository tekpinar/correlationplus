# CorrelationPlus

A Python script to plot and analyze dynamical correlations of proteins.

CorrelationPlus is a script that you can use to plot and analyze 
dynamical correlations for proteins and biological macromolecules. 
These correlations can be dynamical cross-correlations or linear mutual
information. 

This program is more useful if your structure contains multiple
chains. The program will produce an output for overall structure 
and all individual intra-chain correlations, if exist. Moreover,  
the program will give you inter-chain correlations, if you have 
more than one chain. 

The program only requires a pdb file and a correlation data matrix. 
The correlation data has to be in matrix format, where only A(i,j)
values are listed in a square matrix format. 

To run a simple example, go to examples folder and then run:
This will produce plot absolute values of dynamical cross correlations.

##

mapAnalysis -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -s absdcc


Sometimes, we may need to plot difference map of two correlation maps:

Finally, correlationPlus can do centrality analysis for your protein.
It computes degree, closeness, betweenness, current flow closeness, 
current flow betweenness and eigenvector centrality.

python ../src/correlationPlus.py centralityAnalysisApp -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -s absdcc

## Installation

### for users

We recommend to use pip
```pip install correlationPlus```
or if you do not have administration rights
```pip install --user correlationPlus```
If you prefer to use a virtualenv
```
python3.8 -m venv correlationPlus
cd correlationPlus
source bin/activate
pip install correlationPlus
```bash

### for developper

We recommend to use pip and a virtualenv
python3.8 -m venv correlationPlus
cd correlationPlus
source bin/activate
mkdir src
cd src
git clone git@github.com:tekpinar/correlationPlus.git # or 
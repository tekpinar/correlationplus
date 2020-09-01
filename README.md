# CorrelationPlus

A Python API to plot and analyze dynamical correlations of proteins.

CorrelationPlus contains three scripts that you can use to plot and analyze 
dynamical correlations for proteins and biological macromolecules. 
These correlations can be dynamical cross-correlations or linear mutual
information. 

mapAnalysis app plots and produces analysis scripts for 
correlation maps. It can be more useful if your structure contains multiple
chains. The program will produce an output for overall structure 
and all individual intra-chain correlations, if exist. Moreover,  
the program will give you inter-chain correlations, if you have 
more than one chain. The program only requires a pdb file and a 
correlation data matrix. The correlation data has to be in matrix format, 
where only A(i,j) values are listed in a square matrix format. You can 
analyze the correlations with VMD just by loading the tcl files produced by
mapAnalysis module/ 

## A Quick Start with correlationPlus scripts

To run a simple example, go to examples folder and then run:

```bash
correlationPlus mapAnalysis -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -t absndcc
```
This will produce plots of absolute values of dynamical cross correlations.

Sometimes, we may need to plot difference map of two correlation maps. 
For example, you may want to see the differences of linear mutual information 
maps of produced with two different methods, conditions etc.
This can be produced with diffMap app as follows:  

```bash
correlationPlus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi
```

Finally, correlationPlus can do centrality analysis for your protein
via its centralityAnalysis app.

It computes degree, closeness, betweenness, current flow closeness, 
current flow betweenness and eigenvector centrality.

```bash
correlationPlus centralityAnalysis -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -t absndcc
```

## Ipython Interface
For a detailed analysis, script interfaces provided by mapAnalysis, diffMap and 
centralityAnalysis apps may not be sufficient. Therefore, you can use IPython 
to load the functions and do a detailed analysis as follows. 

```
from correlationPlus.mapAnalysis import *
```
 

You can get help for individual functions with

```
help(intraChainCorrelationMaps) 

```
You can different valueFilters, distanceFilters for your analysis. 
Even you can scan a range of values by calling the functions in a 
loop. 



## Installation

### for users

We recommend to use pip
```bash
pip install correlationPlus
```

or if you do not have administration rights
```bash
pip install --user correlationPlus
```

If you prefer to use a virtualenv
```bash
python3.8 -m venv correlationPlus
cd correlationPlus
source bin/activate
pip install correlationPlus
```

### for developers

We recommend to use pip and a virtualenv
```bash
python3.8 -m venv correlationPlus
cd correlationPlus
source bin/activate
mkdir src
cd src
git clone git@github.com:tekpinar/correlationPlus.git # or https://github.com/tekpinar/correlationPlus.git
cd correlationPlus
pip install -e .
```

## Licensing

correplationPlus is developed and released under GNU Lesser GPL Licence. 
Please read to the COPYING and COPYING.LESSER files to know more. 
# CorrelationPlus

A Python package to calculate, visualize and analyze dynamical correlations of proteins.

CorrelationPlus contains four scripts that you can use to calculate, visualize
and analyze dynamical correlations of proteins. 
These correlations can be dynamical cross-correlations, linear mutual
information or generalized correlations. 

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
python3 -m venv correlationPlus
cd correlationPlus
source bin/activate
pip install correlationPlus
```

### for developers

We recommend to use pip and a virtualenv
```bash
python3 -m venv correlationPlus
cd correlationPlus
source bin/activate
mkdir src
cd src
git clone https://github.com/tekpinar/correlationPlus.git # or git@github.com:tekpinar/correlationPlus.git
cd correlationPlus
pip install -e .
```

### from Docker image

Docker images are also available from [Docker Hub](https://hub.docker.com/r/structuraldynamicslab/correlation_plus)

```bash
mkdir shared_dir
cp 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt 6fl9_centeredOrientedAligned2Z.pdb shared_dir
chmod 777 shared_dir
cd shared_dir
docker run -v $PWD:/correlationPlus structuraldynamicslab/correlation_plus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi
```

### from Singularity image

As the docker image is registered in dockerhub you can also use it directly with [Singularity](https://sylabs.io/docs/) 

```bash
singularity run docker://structuraldynamicslab/correlation_plus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi
```
or in 2 steps

```bash
singularity pull correlation_plus.simg docker://structuraldynamicslab/correlation_plus
./correlation_plus.simg -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi
```

Unlike docker you have not to worry about shared directory, your home and /tmp are automatically shared.


## A Quick Start with correlationPlus Scripts

### Calculating dynamical cross-correlations
Download examples folder and go there. 

To calculate dynamical cross-correlations with Gaussian network model:

```bash
correlationPlus calculate -p 6fl9_centeredOrientedAligned2Z.pdb -m GNM -o gnm-ndcc.dat
```

To calculate dynamical cross-correlations with Anisotropic network model:

```bash
correlationPlus calculate -p 6fl9_centeredOrientedAligned2Z.pdb -m ANM -o anm-ndcc.dat
```

### Visualization of correlation maps
To run a simple example of visualization, you can use data and pdb files in the examples folder:

```bash
correlationPlus visualizemap -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -t absndcc
```
This will produce plots of absolute values of dynamical cross correlations.

The visualizemap will produce an output for overall structure 
and all individual intra-chain correlations, if exist. Moreover, the program 
will give you inter-chain correlations, if you have more than one chain. 
The correlation data has to be in matrix format, where only A(i,j) values are 
listed in a square matrix format. LMI matrices produced by g_correlation 
program of Lange and Grubmuller can also be parsed. 

You can analyze the correlations with VMD just by loading the tcl files produced by 
visualizemap script. You can call VMD and go to Extensions->Tk Console menu. 
Write the following command to see the correlations:
```bash
source correlation-interchain-chainsA-B.tcl
```

If you prefer to do the tcl loading in a single command:
```bash
vmd -e correlation-interchain-chainsA-B.tcl
```
Please, beware that the loading can take some time depending on protein size
and number of correlations. 

Additionally, vmd command has to be in your path if you want to do this 
with the command above.
 
Sometimes, we may need to plot difference map of two correlation maps. 
You may want to see the differences of linear mutual information 
maps of produced with two different methods, conditions etc. The correlations
of ligated vs unligated simulations are some common examples.  
The difference maps can be produced with diffMap app as follows:  

```bash
correlationPlus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi
```

### Centrality analysis of the correlation maps
Centrality analysis can be used to deduce active sites, binding sites, 
key mutation sites and allosteric residues. 
correlationPlus can do centrality analysis for your protein
via its centralityAnalysis app.

It can compute degree, closeness, betweenness, current flow closeness, 
current flow betweenness and eigenvector centrality. The following command 
will do all of the above analysis:

```bash
correlationPlus centralityAnalysis -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -t absndcc
```
After the calculation, the centrality values will be inserted into Bfactor
 column of a pdb file. You can load the pdb files with your favorite visualization 
software and color according to Bfactors. If you prefer VMD - as we do-, 
the app will produces tcl files so that you can visualize the key residues with VMD.
The tcl script highlights the residues with the highest 10% of the selected centrality
in VDW representation.

```bash
vmd -e correlation_degree.tcl
```

## Ipython Interface
For a detailed analysis, script interfaces provided by visualizemap, diffMap and 
centralityAnalysis apps may not be sufficient. Therefore, you can use IPython 
to load the functions and do a detailed analysis as follows. 

```
from correlationPlus.visualize import *
```
 

You can get help for individual functions with

```
help(intraChainCorrelationMaps) 

```
You can check different valueFilters, distanceFilters for your analysis. 
Even you can scan a range of values by calling the functions in a 
loop. 


## Licensing

correplationPlus is developed and released under [GNU Lesser GPL Licence](https://www.gnu.org/licenses/lgpl-3.0.en.html). 
Please read to the COPYING and COPYING.LESSER files to know more. 
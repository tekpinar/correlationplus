Quickstart
==========

There are three main scrips of correlationplus package: 
* calculate
* visualize
* analyze

You can find more information about each module as follows:

``correlationplus calculate -h``

``correlationplus visualize -h``

``correlationplus analyze -h``

Here are some examples of the correlationplus commandline interface.
You can find all required files in the examples folder at our `github page <https://github.com/tekpinar/correlationplus>`_

**calculate** module
--------------------
With this module, you can perform dynamical cross-correlation and linear mutual information from
elastic network models and molecular dynamics trajectories. 

To calculate **dynamical cross-correlations** with **Gaussian** network model:

``correlationplus calculate -p 6fl9_centeredOrientedAligned2Z.pdb -m GNM -o gnm-ndcc.dat``

To calculate **dynamical cross-correlations** with **Anisotropic** network model:

``correlationplus calculate -p 6fl9_centeredOrientedAligned2Z.pdb -m ANM -o anm-ndcc.dat``

To calculate **dynamical cross-correlations** from a molecular dynamics trajectory (in dcd, xtc or trr format):

``correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -f 6lu7_dimer_with_N3_protein_sim1_ca_short.trr``

To calculate **linear mutual informations** with **Anisotropic** network model:

``correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi``

To calculate **linear mutual informations** from a molecular dynamics trajectory (in dcd, xtc or trr format):

``correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -f 6lu7_dimer_with_N3_protein_sim1_ca_short.trr -t lmi``

**visualize** module
--------------------
Visualize module plots all 2D correlation maps. It prepares tcl file and the correlated residue pairs can
be visualized with the help of **VMD** program. This interface needs only a pdb file with N residues and
a square matrix of NxN. The correlation data has to be in matrix format, where only A(i,j) values are 
listed in a square matrix format. LMI matrices produced by g_correlation program of Lange and Grubmuller
can also be parsed. 


To run a simple example of visualization, you can use data and pdb files in the examples folder:

``correlationplus visualize -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -t absndcc``

In addition, the command above will produce plots of absolute values of dynamical cross correlations vs interresidue distances.

The visualize app will produce an output for overall structure 
and all individual intra-chain correlations, if exist. Moreover, the program 
will give you inter-chain correlations, if you have more than one chain. 

You can analyze the correlations with `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ just by loading the tcl files produced by 
visualize app. You can call *VMD* and go to *Extensions->Tk Console menu*. 
Write the following command to see the correlations:

``source correlation-interchain-chainsA-B.tcl``

If you prefer to do the tcl loading in a single command:

``vmd -e correlation-interchain-chainsA-B.tcl``

Please, beware that the loading can take some time depending on protein size
and number of correlations. 

Additionally, vmd command has to be in your path if you want to do this 
with the command above.
 
Sometimes, we may need to plot difference map of two correlation maps. 
You may want to see the differences of linear mutual information 
maps of produced with two different methods, conditions etc. The correlations
of ligated vs unligated simulations are some common examples.  
The difference maps can be produced with diffMap app as follows:  

``correlationplus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi``

**analyze** module
------------------
This module can be used to perform centrality analysis of the correlation maps.
Centrality analysis is used to deduce active sites, binding sites, key mutation
sites and allosteric residues. 

The script can compute degree, closeness, betweenness, current flow closeness, 
current flow betweenness and eigenvector centralities. The following command 
will do all of the above analysis:

``correlationplus analyze -i 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt -p 6fl9_centeredOrientedAligned2Z.pdb -t absndcc``

After the calculation, the centrality values will be inserted into **Bfactor** 
column of a pdb file. You can load the pdb files with your favorite visualization 
software and color according to **Bfactors**. If you prefer *VMD* - as we do-, 
the app will produce tcl files so that you can visualize the key residues with **VMD**.
The tcl script highlights the residues with the highest 10% of the selected centrality
in VDW representation.

``vmd -e correlation_degree.tcl``

Ipython Interface
-----------------
For a detailed analysis, script interfaces provided by calculate, visualize, analyze and 
diffMap apps may not be sufficient. Therefore, you can use IPython 
to load the functions and do a detailed analysis as follows. 

``from correlationplus.visualize import *``
 
You can get help for each function with

``help(intraChainCorrelationMaps)``

You can check different valueFilters, distanceFilters for your analysis. 
Also, you can scan a range of values by calling the functions in a 
loop. 

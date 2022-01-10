Quickstart
==========

There are four main scrips of correlationplus package:

* calculate
* visualize
* analyze
* paths

You can find more information about each script as follows::

    correlationplus calculate -h

    correlationplus visualize -h

    correlationplus analyze -h
    
    correlationplus paths -h

Here are some examples of the correlationplus commandline interface.
You can find all required files in the examples folder at our `github page <https://github.com/tekpinar/correlationplus>`_

Correlation Data Types
----------------------
Correlationplus can handle the following correlation/coupling data types:

* dcc: Dynamical cross-correlations in full matrix format.
* ndcc: Normalized dynamical cross-correlations in full matrix format.
* absndcc: Absolute values normalized dynamical cross-correlations in full matrix format
* omegacc: Normalized Pearson correlations of backbone dihedral angle omega in full matrix format.
* phicc: Normalized Pearson correlations of backbone dihedral angle phi in full matrix format.
* psicc: Normalized Pearson correlations of backbone dihedral angle psi in full matrix format.
* lmi: Linear mutual information in full matrix format or output of g_correlation program.
* nlmi: Normalized linear mutual information in full matrix format or output of g_correlation program. 
* coeviz: After removing the header lines, the data is in matrix format. 
* evcouplings: The sequence coupling csv files obtained from https://evcouplings.org/ can be parsed directly. 
* generic: If you have some coupling data (from dynamics, sequences or any other data) in full matrix format, select this option. 

You can control data types in correlationplus script with '-t' option. For example, '-t generic' will tell the script that it is 
a generic data, etc.


**calculate** script
--------------------
With this module, you can calculate dynamical cross-correlation and linear mutual information from
elastic network models and molecular dynamics trajectories. 

The only input file needed is a PDB file. The file can contain single or multiple chains. In all of 
the computations via script interfaces, only Calpha atoms are selected and used.    

To calculate **normalized dynamical cross-correlations** with **Gaussian** network model::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb -m GNM -o ndcc-6lu7-gnm.dat

To calculate **normalized dynamical cross-correlations** with **Anisotropic** network model::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb -m ANM -o ndcc-6lu7-anm.dat

To calculate **normalized dynamical cross-correlations** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			                      -o ndcc-6lu7-md.dat

To calculate **normalized linear mutual informations** with **Anisotropic** network model::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb -t nlmi -o nlmi-6lu7-anm.dat

To calculate **normalized linear mutual informations** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr -t nlmi\
			                      -o nlmi-6lu7-md.dat
To calculate **normalized Pearson cross-correlations of backbone omega dihedral angles** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			    -t omegacc -o omegacc-6lu7-md.dat

To calculate **normalized Pearson cross-correlations of backbone phi dihedral angles** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			    -t phicc -o phicc-6lu7-md.dat

To calculate **normalized Pearson cross-correlations of backbone psi dihedral angles** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			    -t psicc -o psicc-6lu7-md.dat

Sometimes, there are not dihedral angles for some residues at the beginning/end of the chains (See https://userguide.mdanalysis.org/1.1.1/examples/analysis/structure/dihedrals.html for details). If there are some missing atoms, you may not also be able to calculate some dihedral angles. To avoid these problems, we fill these missing values with 1.0 to maintain a one-to-one correspondence with the number of residues. Due to this reason, pay attention to the highly correlated values that may be due to this artificial filling. 

**visualize** script
--------------------
Visualize module plots all 2D correlation maps. It prepares tcl and pml files so that the correlated residue pairs can
be visualized with the help of **VMD** and **PyMOL** programs. This interface needs only a pdb file with N residues and
a square matrix of NxN. The correlation data has to be in matrix format, where only A(i,j) values are 
listed in a square matrix format. LMI matrices produced by g_correlation program of Lange and Grubmuller
can also be parsed. 


To run a simple example of visualization, you can use the data and pdb files in the examples folder::

  correlationplus visualize -i ndcc-6lu7-anm.dat -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t absndcc -v 0.75

In addition, the command above will produce plots of absolute values of dynamical cross correlations vs interresidue distances.
This information can be quite useful if you are particularly looking for long-distance interactions. 

The visualize app will produce an output for overall structure and all individual intra-chain correlations, if exist. 
Moreover, the program will give you inter-chain correlations, if you have more than one chain. 

You can analyze the correlations with `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ just by loading the tcl files produced by 
visualize app.  To reduce the clutter, the command above will only dump the correlations greater than 0.75 to your tcl or pml file.
If you would like to visualize an interval, you can specify the maximal value as well with '-x ' parameter.

You can call VMD and go to *Extensions->Tk Console menu*. 
Write the following command to see the correlations::

  source correlation-interchain-chainsA-B.tcl

If you prefer to do the tcl loading in a single command::

  vmd -e correlation-interchain-chainsA-B.tcl

Please, beware that the loading can take some time depending on protein size,
number of correlations and the min-max correlation limits that you imposed. 

Additionally, vmd command has to be in your path if you want to do this 
with the command above.

If you would like to use PyMOL, the following command will be sufficient::
  
  pymol correlation-interchain-chainsA-B.pml



Sometimes, we may need to plot difference map of two correlation maps. 
You may want to see the differences of linear mutual information 
maps of produced with two different methods, conditions etc. The correlations
of ligated vs unligated simulations are some common examples.  
The difference maps can be produced with diffMap app as follows::

  correlationplus diffMap -i 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
                          -j 6lu7_dimer_no_ligand_protein_sim1-lmi.dat\
			  -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

**analyze** script
------------------
This module can be used to perform centrality analysis of the correlation maps.
Centrality analysis is used to deduce active sites, binding sites, key mutation
sites and allosteric residues. 

The script can compute degree, closeness, betweenness, current flow closeness, 
current flow betweenness, eigenvector centralities and major communities. The following
command will do all of the above analysis::

  correlationplus analyze -i 6lu7_dimer_with_N3_protein_sim1-lmi.dat\
                          -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

If you would like to calculate only a certain centrality like betweenness::

  correlationplus analyze -i 6lu7_dimer_with_N3_protein_sim1-lmi.dat\
                          -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb
			  -c betweenness -t lmi

After the calculation, the centrality values will be inserted into **Bfactor** 
column of a pdb file. You can load the pdb files with your favorite visualization 
software and color according to **Bfactors**. If you prefer **VMD** - as we do-, 
the app will produce tcl files so that you can visualize the key residues with **VMD**.
The tcl script highlights the residues with the highest 10% of the selected centrality
in VDW representation.::

  vmd -e correlation_degree.tcl

With PyMol::
  
  pymol correlation_degree.pml

**paths** script
------------------
To calculate suboptimal paths between two active site residues in chain A and chain B of 
SARS-CoV2 main protease::

    correlationplus paths -i ndcc-6lu7-anm.dat\
              		  -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb\
              		  -b A41 -e B41
   
This command will only produce the optimal path and print out the path length. If you would like
to calculate suboptimal paths as well, you can append -n argument. Here is the example command to 
calculate 10 paths between residue 41 of chain A and residue 41 of chain B::

    correlationplus paths -i ndcc-6lu7-anm.dat\
              		  -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb\
              		  -b A41 -e B41 -n 10



Ipython Interface
-----------------
For a detailed analysis, script interfaces provided by calculate, visualize, analyze, paths and 
diffMap scripts may not be sufficient. Therefore, you can use IPython 
to load the modules and do a detailed analysis as follows. 

``from correlationplus.calculate import *``

``from correlationplus.visualize import *``
 
You can get help for each function with

``help(intraChainCorrelationMaps)``

You can check different valueFilters, distanceFilters for your analysis. 
Also, you can scan a range of values by calling the functions in a 
loop. 

There is a minor but important difference between the scripts and the modules for centrality and path 
analyses. If you want to use the module for centrality analysis:

``from correlationplus.centralityAnalysis import *``

Please notice that the name of the script was 'analyze' but the name of the module is 'centralityAnalysis'. 

Similarly, the name of the path analysis script is 'paths' while the name of the module is 'pathAnalysis'. 
Therefore, you have to call path analysis module interactively as follows:

``from correlationplus.pathAnalysis import *``






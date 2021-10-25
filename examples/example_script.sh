#!/bin/bash
#Purpose     : To show basic usage of correlationplus Python package on a Unix system.#              
#              The script has been tested  with 0.2.0 version of correlationplus on
#              CentOS 7 and MacOS Catalina.
#              Please visit https://correlationplus.readthedocs.io/en/latest/quickstart.html
#              for the most uptodate usage information and more examples. 
#Prerequisite: You must have installed correlationplus. Otherwise, try 'pip install correlationplus'.
#Author      : Mustafa Tekpinar
#Email       : tekpinar@buffalo.edu
#Date        : August 9, 2021
#Licence     : LGPL v3

#We recommend to run it in a completely empty folder.
#The script will try to download a monomer of apo SARS-CoV-2 main protease pdb file 6y2e.pdb
#with curl. If you have wget, try to replace curl with wget.

#Step 0: Download 6y2e.pdb with curl or wget.
#If you want, you can download your own pdb as well.
#In this case, do not forget to comment wget and curl lines.
#In addition, correct the protein variable in the following line.
protein="6y2e"


#wget https://files.rcsb.org/download/${protein}.pdb
curl https://files.rcsb.org/download/${protein}.pdb -o ${protein}.pdb

if test -f "${protein}.pdb"
then
    echo "${protein}.pdb was downloaded successfully."

    #Step 1: Calculate normalized LMI matrix from CG-ANM
    echo "Calculating normalized LMI matrix from CG-ANM."
    correlationplus calculate -p ${protein}.pdb -t nlmi -o nLMI.dat

    #Step 2: Visualize the normalized LMI matrix
    echo "Visualizing normalized LMI matrix"
    echo "Only high correlations (>=0.75) will be projected onto the protein!"
    correlationplus visualize -p ${protein}.pdb -i nLMI.dat -t nlmi -o nlmi -v 0.75

    #Step 3: Network analysis of the normalized LMI matrix
    echo "Analyzing dynamical network from the normalized LMI matrix"
    correlationplus analyze -p ${protein}.pdb -i nLMI.dat -t lmi -o network

    #Voila! You have many pml and tcl files.
    #You can open pml files as
    #pymol network_eigenvector.pml

    #If you are a VMD person
    #vmd -e network_eigenvector.tcl

    #Please note that pymol or vmd has to be in your path to run these commands.
    #Any other problem?
    #Read the documentation at https://correlationplus.readthedocs.io/en/latest/quickstart.html
    #Or email me at tekpinar@buffalo.edu!
    
else    
        echo "ERROR: ${protein}.pdb file is not in the directory!"
	echo "wget or curl failed to download the protein!"
fi

##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   correlation maps of proteins.                             #
# Authors: Mustafa Tekpinar                                                   #
# Copyright (C) Mustafa Tekpinar, 2017-2018                                   #
# Copyright (C) CNRS-UMR3528, 2019                                            #
# Copyright (C) Institut Pasteur Paris, 2020-2021                             #
#                                                                             #
# This file is part of correlationplus.                                       #
#                                                                             #
# correlationplus is free software: you can redistribute it and/or modify     #
# it under the terms of the GNU Lesser General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# correlationplus is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU LESSER General Public License for more details.                         #
#                                                                             #
# You should have received a copy of the GNU Lesser General Public License    #
# along with correlationplus.  If not, see <https://www.gnu.org/licenses/>.   #
###############################################################################


import sys
import getopt

import numpy as np
from prody import parsePDB
from prody import buildDistMatrix
from correlationplus.visualize import convertLMIdata2Matrix
from correlationplus.visualize import parseEVcouplingsScores
from correlationplus.visualize import parseSparseCorrData
from correlationplus.visualize import parseElasticityGraph
from correlationplus.centralityAnalysis import centralityAnalysis
from correlationplus.centralityAnalysis import buildDynamicsNetwork
from correlationplus.centralityAnalysis import buildSequenceNetwork
# from .visualizemap import handle_arguments_visualizemapApp


def usage_centralityAnalysisApp():
    """
    Show how to use this program!
    """
    print("""
Example usage:
correlationplus analyze -i ndcc-6lu7-anm.dat -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb

Arguments: -i: A file containing correlations in matrix format. (Mandatory)

           -p: PDB file of the protein. (Mandatory)
           
           -t: Type of the matrix. It can be ndcc, lmi, absndcc (absolute values of ndcc)
               or eg (elasticity graph).
               In addition, coeviz and evcouplings are also some options to analyze sequence
               correlations. 
               If your data any other coupling data in full matrix format, you can select generic
               as your data type. 
               Default value is ndcc (Optional)

           -o: This will be your output file. Output figures are in png format. 
               (Optional)

           -c: Type of the centrality that you want to calculate. Default is 'all'.
               (Optional). If you want only a particular centrality, you can select 
               one of the following options: 
                    * degree
                    * betweenness
                    * closeness
                    * current_flow_betweenness
                    * current_flow_closeness
                    * eigenvector
                    * community
               
           -v: Value filter. The values lower than this value in the map will be 
               considered as zero. Default is 0.3. (Optional)

           -d: Distance filter. The residues with distances higher than this value 
               will be considered as zero. Default is 7.0 Angstrom. (Optional)
""")


def handle_arguments_centralityAnalysisApp():
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None
    centrality_type = None
    value_cutoff = None
    distance_cutoff = None

    try:
        opts, args = getopt.getopt(sys.argv[2:], "hi:o:t:p:c:v:d:", ["help", "inp=", "out=", "type=", "pdb=", "centrality=", "value=", "distance="])
    except getopt.GetoptError:
        usage_centralityAnalysisApp()
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            usage_centralityAnalysisApp()
            sys.exit(-1)
        elif opt in ("-i", "--inp"):
            inp_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-t", "--type"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file = arg
        elif opt in ("-c", "--centrality"):
            centrality_type = arg
        elif opt in ("-v", "--value"):
            value_cutoff = arg
        elif opt in ("-d", "--distance"):
            distance_cutoff = arg
        else:
            assert False, usage_centralityAnalysisApp()

    # Input data matrix and PDB file are mandatory!
    if inp_file is None or pdb_file is None:
        print("@> ERROR: A PDB file and a correlation matrix are mandatory!")
        usage_centralityAnalysisApp()
        sys.exit(-1)

    # Assign a default name if the user forgets the output file prefix.
    if out_file is None:
        out_file = "correlation"

    # The user may prefer not to submit a title for the output.
    if sel_type is None:
        sel_type = "ndcc"

    if centrality_type is None:
        centrality_type = "all"

    if value_cutoff is None:
        value_cutoff = 0.3

    if distance_cutoff is None:
        distance_cutoff = 7.0

    return inp_file, out_file, sel_type, pdb_file, centrality_type, value_cutoff, distance_cutoff


def centralityAnalysisApp():
    inp_file, out_file, sel_type, pdb_file, centrality_type, value_cutoff,\
            distance_cutoff = handle_arguments_centralityAnalysisApp()

    print(f"""
@> Running 'analyze' app

@> Input file     : {inp_file}
@> PDB file       : {pdb_file}
@> Data type      : {sel_type}
@> Output         : {out_file}
@> Centrality     : {centrality_type}
@> Value filter   : {value_cutoff}
@> Distance filter: {distance_cutoff}""")

    ##########################################################################
    # Read PDB file
    # TODO: This is the only place where I use Prody.
    # Maybe, I can replace it with a library that only parses
    # PDB files. Prody does a lot more!
    selectedAtoms = parsePDB(pdb_file, subset='ca')
    valueFilter = float(value_cutoff)
    distanceFilter = float(distance_cutoff)
    distanceMatrix = buildDistMatrix(selectedAtoms)

    ##########################################################################
    # Read data file and assign to a numpy array
    if sel_type.lower() == "ndcc":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()
 
        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True, 
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = np.loadtxt(inp_file, dtype=float)
    elif sel_type.lower() == "absndcc":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()
 
        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = np.absolute(parseSparseCorrData(inp_file, selectedAtoms, \
                                                        Ctype=True, 
                                                        symmetric=True,
                                                        writeAllOutput=False))
        else:
            ccMatrix = np.absolute(np.loadtxt(inp_file, dtype=float))
    elif sel_type.lower()== "lmi":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()
 
        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True, 
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=False)
    elif sel_type.lower() == "coeviz":
        ccMatrix = np.loadtxt(inp_file, dtype=float) 
    elif sel_type.lower() == "evcouplings":
        ccMatrix = parseEVcouplingsScores(inp_file, selectedAtoms, False)
    elif sel_type.lower() == "generic":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()
 
        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True, 
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = np.loadtxt(inp_file, dtype=float)
    elif sel_type.lower() == "eg":
        # The data type is elasticity graph
        ccMatrix = parseElasticityGraph(inp_file, selectedAtoms, \
                                            writeAllOutput=False)
    else:
        print("@> ERROR: Unknown data type: Type can only be ndcc, absndcc, lmi,\n")
        print("@>        coeviz or evcouplings. If you have your data in full \n")
        print("@>        matrix format and your data type is none of the options\n")
        print("@>        mentionned, you can set data type 'generic'.\n")
        sys.exit(-1)

    if ((sel_type.lower() == "evcouplings") or \
        (sel_type.lower() == "generic")  or \
        (sel_type.lower() == "eg")):
        network = buildSequenceNetwork(ccMatrix, distanceMatrix, \
                                    valueFilter, distanceFilter,\
                                    selectedAtoms)
    else:
        network = buildDynamicsNetwork(ccMatrix, distanceMatrix, \
                                    valueFilter, distanceFilter,\
                                    selectedAtoms)

    if centrality_type == "all":
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, "degree",
                           selectedAtoms)
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, "betweenness",
                           selectedAtoms)
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, "closeness",
                           selectedAtoms)
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, "current_flow_betweenness",
                           selectedAtoms)
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, "current_flow_closeness",
                           selectedAtoms)
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, "eigenvector",
                           selectedAtoms)
        # Community analysis is time consuming. Therefore, it will not be called by default.
        # centralityAnalysis(ccMatrix, valueFilter, distanceFilter, out_file, "community",
        #                    selectedAtoms)
    else:
        centralityAnalysis(network, valueFilter, distanceFilter, out_file, centrality_type,
                           selectedAtoms)

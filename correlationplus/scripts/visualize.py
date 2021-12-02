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
import os

#TODO: Move to argparse instead of getopt
import getopt
from collections import Counter

import numpy as np
from prody import parsePDB

from correlationplus.visualize import *

def usage_visualizemapApp():
    """
    Show how to use this program!
    """
    print("""
Example usage:
correlationplus visualize -i ndcc-6lu7-anm.dat -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb

Arguments: -i: A file containing correlations in matrix format. (Mandatory)

           -p: PDB file of the protein. (Mandatory)
           
           -t: Type of the matrix. It can be dcc, ndcc, lmi, nlmi (normalized lmi), 
               absndcc (absolute values of ndcc) or eg (elasticity graph).
               In addition, coeviz and evcouplings are also some options to analyze sequence
               correlations. 
               If your data is any other coupling data in full matrix format, you can select 
               your data type as 'generic'. 
               Default value is ndcc (Optional)

           -v: Minimal correlation value. Any value equal or greater than this 
               will be projected onto protein structure. Default is minimal value of the map. (Optional)
            
           -x: Maximal correlation value. Any value equal or lower than this 
               will be projected onto protein structure. Default is maximal value of the map. (Optional)

           -d: If the distance between two C_alpha atoms is bigger 
               than the specified distance, it will be projected onto protein structure. 
               Default is 0.0 Angstroms. (Optional)

           -D: If the distance between two C_alpha atoms is smaller 
               than the specified distance, it will be projected onto protein structure. 
               Default is 9999.9 Angstroms. (Optional)

           -r: Cylinder radius scaling coefficient to multiply with the correlation quantity.
               It can be used to improve tcl and pml outputs to view the interaction 
               strengths properly. Recommended values are between 0.0 and 2.0. (Optional)

           -o: This will be your output file. Output figures are in png format. (Optional)
""")


def handle_arguments_visualizemapApp():
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None
    vmin_fltr = None
    vmax_fltr = None
    dmin_fltr = None
    dmax_fltr = None
    cyl_rad = None
    try:
        opts, args = getopt.getopt(sys.argv[2:], "hi:o:t:p:v:x:d:D:r:", \
            ["help", "inp=", "out=", "type=", "pdb=", "vmin=", "vmax=", \
            "dmin=", "dmax=", "radius="])
    except getopt.GetoptError:
        usage_visualizemapApp()
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            usage_visualizemapApp()
            sys.exit(-1)
        elif opt in ("-i", "--inp"):
            inp_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-t", "--type"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file = arg
        elif opt in ("-v", "--vmin"):
            vmin_fltr = arg
        elif opt in ("-x", "--vmax"):
            vmax_fltr = arg
        elif opt in ("-d", "--dmin"):
            dmin_fltr = arg
        elif opt in ("-D", "--dmax"):
            dmax_fltr = arg
        elif opt in ("-r", "--radius"):
            cyl_rad = arg
        else:
            assert False, usage_visualizemapApp()

    # Input data matrix and PDB file are mandatory!
    if inp_file is None or pdb_file is None:
        usage_visualizemapApp()
        sys.exit(-1)

    if not os.path.exists(inp_file):
        print("@> ERROR: Could not find " + inp_file)
        sys.exit(-1)
    if not os.path.exists(pdb_file):
        print("@> ERROR: Could not find " + pdb_file)
        sys.exit(-1)
    # Assign a default name if the user forgets the output file prefix.
    if out_file is None:
        out_file = "correlation"

    # The user may prefer not to submit a title for the output.
    if sel_type is None:
        sel_type = "ndcc"

    # if vmin_fltr is None:
    #     vmin_fltr = 0.75

    if dmin_fltr is None:
        dmin_fltr = 0.0
    
    if dmax_fltr is None:
        dmax_fltr = 9999.9

    return inp_file, out_file, sel_type, pdb_file, \
            vmin_fltr, vmax_fltr, \
            dmin_fltr, dmax_fltr, cyl_rad


def visualizemapApp():
    inp_file, out_file, sel_type, pdb_file, \
    vmin_fltr, vmax_fltr, \
    dmin_fltr, dmax_fltr, cyl_rad = handle_arguments_visualizemapApp()
    print(f"""
@> Running 'visualize' app:
    
@> Input file       : {inp_file}
@> PDB file         : {pdb_file}
@> Data type        : {sel_type}
@> Min. dist. filter: {dmin_fltr}
@> Max. dist. filter: {dmax_fltr}
@> Output           : {out_file}""")


    ##########################################################################
    # Read PDB file
    protein= parsePDB(pdb_file)
    selectedAtoms = protein.select('protein and name CA')
    print("@> Selected only "+str(selectedAtoms.numAtoms())+" atoms!\n")
    ##########################################################################
    minColorBarLimit = 0.0
    maxColorBarLimit = 1.0
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
        # Check the data range in the matrix.
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        if minCorrelationValue < 0.0:
        # Assume that it is an nDCC file
            minColorBarLimit = -1.0
        if maxCorrelationValue > 1.00001:
            print("This correlation map is not normalized!")
            # TODO: At this point, one can ask the user if s/he wants to normalize it!
            sys.exit(-1)
        else:
            maxColorBarLimit = 1.0
    elif sel_type.lower() == "dcc":
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
        # Check the data range in the matrix.
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
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
        minColorBarLimit = 0.0
        maxColorBarLimit = 1.0
    elif sel_type.lower() == "lmi":
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
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
        #minColorBarLimit = 0.0

    elif sel_type.lower() == "nlmi":
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
        #minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = 0.0

        # Ideally, it is supposed to be 1 but I used 1.00001 to avoid  
        # rounding problems
        if maxCorrelationValue > 1.00001:
            print("This LMI map is not normalized!")
            # TODO: At this point, one can ask the user if s/he wants to normalize it!
            sys.exit(-1)
        else:
            maxColorBarLimit = 1.0

    elif sel_type.lower() == "coeviz":
        ccMatrix = np.loadtxt(inp_file, dtype=float)
        minColorBarLimit = 0.0
        maxColorBarLimit = 1.0
    
    elif sel_type.lower() == "evcouplings":
        ccMatrix = parseEVcouplingsScores(inp_file, selectedAtoms, False)
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
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

        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
    elif sel_type.lower() == "eg":
        # The data type is elasticity graph
        ccMatrix = parseElasticityGraph(inp_file, selectedAtoms, \
                                            writeAllOutput=False)

        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
    else:
        print("@> ERROR: Unknown data type: Type can only be ndcc, absndcc, lmi,\n")
        print("@>        coeviz or evcouplings. If you have your data in full \n")
        print("@>        matrix format and your data type is none of the options\n")
        print("@>        mentionned, you can set data type 'generic'.\n")
        sys.exit(-1)

    # Set vmin_fltr and vmax_fltr
    if (vmin_fltr == None):
        vmin_fltr = minColorBarLimit
    if (vmax_fltr == None):
        vmax_fltr = maxColorBarLimit
    
    print(f"""@> Min. value filter: {vmin_fltr}""")
    print(f"""@> Max. value filter: {vmax_fltr}""")
    
    ##########################################################################
    # Call overall correlation calculation

    overallCorrelationMap(ccMatrix, minColorBarLimit, maxColorBarLimit,
                          out_file, " ", selectedAtoms)

    plotDistributions = True
    VMDcylinderRadiusScale = 0.5
    PMLcylinderRadiusScale = 0.3

    if (cyl_rad == None):
        if sel_type.lower() == "evcouplings":
            VMDcylinderRadiusScale = 0.02
            PMLcylinderRadiusScale = 0.02
        else:
            VMDcylinderRadiusScale = 0.5
            PMLcylinderRadiusScale = 0.3
        print(f"""@> VMD Cylinder radius: {VMDcylinderRadiusScale}""")
        print(f"""@> PyMol Cylinder radius: {PMLcylinderRadiusScale}""")

    else:
        VMDcylinderRadiusScale = float(cyl_rad)
        PMLcylinderRadiusScale = float(cyl_rad)
        print(f"""@> Cylinder radius: {cyl_rad}""")

    if plotDistributions:
        if sel_type.lower() == "ndcc":
            distanceDistribution(ccMatrix, out_file, "nDCC", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
        elif sel_type.lower() == "dcc":
            distanceDistribution(ccMatrix, out_file, "DCC", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)

        elif sel_type.lower() == "absndcc":
            distanceDistribution(ccMatrix, out_file, "Abs(nDCC)",
                                 selectedAtoms, absoluteValues=True, writeAllOutput=False)

        elif sel_type.lower() == "lmi":
            distanceDistribution(ccMatrix, out_file, "LMI", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)
        elif sel_type.lower() == "nlmi":
            distanceDistribution(ccMatrix, out_file, "nLMI", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)        
        elif sel_type.lower() == "coeviz":
            distanceDistribution(ccMatrix, out_file, "CoeViz", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)

        elif sel_type.lower() == "evcouplings":
            distanceDistribution(ccMatrix, out_file, "EVcoupling Score", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)

        elif sel_type.lower() == "generic":
            distanceDistribution(ccMatrix, out_file, "Correlation", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
        elif sel_type.lower() == "eg":
            distanceDistribution(ccMatrix, out_file, "Force Constants", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
        else:
            print("Warning: Unknows correlation data.\n")
            print("         Correlations can be dcc, ndcc, absndcc, lmi,\n")
            print("         nlmi, coeviz or evcouplings!\n")

    ##########################################################################
    # Check number of chains. If there are multiple chains, plot inter and
    # intra chain correlations
    chains = Counter(selectedAtoms.getChids()).keys()
    saveMatrix = False
    plotChains = True
    if len(chains) > 1 and plotChains:
        intraChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,
                                  out_file, " ", selectedAtoms, saveMatrix)
        interChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,
                                  out_file, " ", selectedAtoms, saveMatrix)

    # Here, we can filter some correlation values closer than a distance.
    # Typically, it is supposed to filter out the correlation within the
    # same secondary structure etc.
    filterByDistance = True
    if filterByDistance:
        disMinValue = float(dmin_fltr)
        disMaxValue = float(dmax_fltr)
        ccMatrix = filterCorrelationMapByDistance(ccMatrix, out_file, " ",
                                                  selectedAtoms, 
                                                  disMinValue, disMaxValue,
                                                  absoluteValues=False,
                                                  writeAllOutput=False)

    # Overall projection
    projectCorrelationsOntoProteinVMD(pdb_file, ccMatrix, out_file,
                                      selectedAtoms, 
                                      vminFilter=float(vmin_fltr),
                                      vmaxFilter=float(vmax_fltr),
                                      cylinderRadiusScaler=VMDcylinderRadiusScale,
                                      absoluteValues=True, writeAllOutput=True)


    projectCorrelationsOntoProteinPyMol(pdb_file, ccMatrix, out_file,
                                      selectedAtoms, 
                                      vminFilter=float(vmin_fltr),
                                      vmaxFilter=float(vmax_fltr),
                                      cylinderRadiusScaler=PMLcylinderRadiusScale,
                                      absoluteValues=True, writeAllOutput=True)

    print("\n@> Program finished successfully!\n")


if __name__ == "__main__":
    visualizemapApp()

###############################################################################
# correlationplus - Python module to plot dynamical correlations maps         #
#                   for proteins.                                             #
# Authors: Mustafa Tekpinar                                                   #
# Copyright (C) Mustafa Tekpinar 2017-2018                                    #
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

from correlationplus.visualize import overallCorrelationMap, convertLMIdata2Matrix, distanceDistribution
from correlationplus.visualize import intraChainCorrelationMaps, interChainCorrelationMaps
from correlationplus.visualize import filterCorrelationMapByDistance, projectCorrelationsOntoProteinVMD
from correlationplus.visualize import projectCorrelationsOntoProteinPyMol
from correlationplus.visualize import parseEVcouplingsScores


def usage_visualizemapApp():
    """
    Show how to use this program!
    """
    print("""
Example usage:
correlationplus visualize -i 4z90-cross-correlations.txt -p 4z90.pdb

Arguments: -i: A file containing correlations in matrix format. (Mandatory)

           -p: PDB file of the protein. (Mandatory)
           
           -t: Type of the matrix. It can be ndcc, lmi or absndcc (absolute values of ndcc).
               In addition, coeviz and evcouplings are also some options to analyze sequence
               correlations. 
               Default value is ndcc (Optional)

           -v: Minimal correlation value. Any value equal or greater than this 
               will be projected onto protein structure. Default is 0.75. (Optional)

           -d: If the minimal distance between two C_alpha atoms is bigger 
               than the specified distance, it will be projected onto protein structure. 
               Default is 0. (Optional)
               
           -o: This will be your output file. Output figures are in png format. (Optional)
""")


def handle_arguments_visualizemapApp():
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None
    vmin_fltr = None
    dis_fltr = None
    try:
        opts, args = getopt.getopt(sys.argv[2:], "hi:o:t:p:v:d:", ["help", "inp=", "out=", "type=", "pdb=", "vmin=", "dis="])
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
        elif opt in ("-d", "--dis"):
            dis_fltr = arg
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

    if vmin_fltr is None:
        vmin_fltr = 0.75

    if dis_fltr is None:
        dis_fltr = 0.0
    

    return inp_file, out_file, sel_type, pdb_file, vmin_fltr, dis_fltr


def visualizemapApp():
    inp_file, out_file, sel_type, pdb_file, vmin_fltr, dis_fltr = \
        handle_arguments_visualizemapApp()
    print(f"""
@> Running 'visualize' app:
    
@> Input file       : {inp_file}
@> PDB file         : {pdb_file}
@> Data type        : {sel_type}
@> Min. value filter: {vmin_fltr}
@> Distance filter  : {dis_fltr}
@> Output           : {out_file}""")


    ##########################################################################
    # Read PDB file
    # TODO: This is the only place where I use Prody.
    # Maybe, I can replace it with a library that only parses
    # PDB files. Prody does a lot more!
    selectedAtoms = parsePDB(pdb_file, subset='ca')

    ##########################################################################
    # Read data file and assign to a numpy array
    if sel_type.lower() == "ndcc":
        ccMatrix = np.loadtxt(inp_file, dtype=float)
        # Check the data range in the matrix.
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        if minCorrelationValue < 0.0:
        # Assume that it is an nDCC file
            minColorBarLimit = -1.0
        if maxCorrelationValue > 1.0:
            print("This correlation map is not normalized!")
            # TODO: At this point, one can ask the user if s/he wants to normalize it!
            sys.exit(-1)
        else:
            maxColorBarLimit = 1.0
    elif sel_type.lower() == "absndcc":
        ccMatrix = np.absolute(np.loadtxt(inp_file, dtype=float))
        minColorBarLimit = 0.0
        maxColorBarLimit = 1.0
    elif sel_type.lower() == "lmi":
        ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=True)
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = 0.0
        if maxCorrelationValue > 1.0:
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
    else:
        print("Unknown matrix data type: The type can only be ndcc, absndcc or lmi!\n")
        sys.exit(-1)

    ##########################################################################
    # Call overall correlation calculation

    overallCorrelationMap(ccMatrix, minColorBarLimit, maxColorBarLimit,
                          out_file, " ", selectedAtoms)

    plotDistributions = True
    VMDcylinderRadiusScale = 0.5
    PMLcylinderRadiusScale = 0.3
    if plotDistributions:
        if sel_type.lower() == "ndcc":
            distanceDistribution(ccMatrix, out_file, "nDCC", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)

        elif sel_type.lower() == "absndcc":
            distanceDistribution(ccMatrix, out_file, "Abs(nDCC)",
                                 selectedAtoms, absoluteValues=True, writeAllOutput=False)

        elif sel_type.lower() == "lmi":
            distanceDistribution(ccMatrix, out_file, "LMI", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)
        
        elif sel_type.lower() == "coeviz":
            distanceDistribution(ccMatrix, out_file, "CoeViz", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)

        elif sel_type.lower() == "evcouplings":
            distanceDistribution(ccMatrix, out_file, "EVcoupling Score", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
            VMDcylinderRadiusScale = 0.01
            PMLcylinderRadiusScale = 0.01
        else:
            print("Warning: Unknows correlation data.\n")
            print("         Correlations can be ndcc, absndcc, lmi,\n")
            print("         coeviz or evcouplings!\n")

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
        distanceValue = float(dis_fltr)
        ccMatrix = filterCorrelationMapByDistance(ccMatrix, out_file, " ",
                                                  selectedAtoms, distanceValue,
                                                  absoluteValues=False,
                                                  writeAllOutput=False)

    # Overall projection
    projectCorrelationsOntoProteinVMD(pdb_file, ccMatrix, out_file,
                                      selectedAtoms, valueFilter=float(vmin_fltr),
                                      cylinderRadiusScaler=VMDcylinderRadiusScale,
                                      absoluteValues=True, writeAllOutput=True)


    projectCorrelationsOntoProteinPyMol(pdb_file, ccMatrix, out_file,
                                      selectedAtoms, valueFilter=float(vmin_fltr),
                                      cylinderRadiusScaler=PMLcylinderRadiusScale,
                                      absoluteValues=True, writeAllOutput=True)

    print("\n@> Program finished successfully!\n")


if __name__ == "__main__":
    visualizemapApp()

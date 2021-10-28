###############################################################################
# correlationplus - Python module to plot dynamical correlations maps         #
#                   for proteins.                                             #
# Authors: Mustafa Tekpinar                                                   #
# Copyright Mustafa Tekpinar 2017-2018                                        #
# Copyright CNRS-UMR3528, 2019                                                #
# Copyright Institut Pasteur Paris, 2020-2021                                 #
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
from collections import Counter

import numpy as np
from prody import parsePDB

from correlationplus.visualize import convertLMIdata2Matrix, intraChainCorrelationMaps, interChainCorrelationMaps
#from correlationplus.diffMap import overallUniformDifferenceMap, overallNonUniformDifferenceMap
from correlationplus.visualize import overallUniformDifferenceMap, overallNonUniformDifferenceMap


def usage_diffMaps():
    """
    Show how to use this program!
    """
    print("""
Example minimal usage:
correlationplus diffMap -i 4z90-cross-correlations.txt -j 4z91-cross-correlations.txt -p 4z90.pdb

Arguments: -i: The first file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)
           -j: The second file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)
           -p: PDB file of the protein. (Mandatory)")
           -q: A second PDB file for the other conformation if residues numbers are not same in two conformations. (Optional)
           -t: It can be ndcc, nlmi or absndcc (absolute values of ndcc). Default value is ndcc (Optional)
           -o: This will be your output file prefix. Output figures are in png format. (Optional)

""")


def handle_arguments_diffMaps():
    inp_file1 = None
    inp_file2 = None
    pdb_file1 = None
    pdb_file2 = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[2:], "hi:j:o:t:p:q:", ["inp1=", "inp2=", "out=", "type=", "pdb=", "pdb2="])
    except getopt.GetoptError:
        usage_diffMaps()
    for opt, arg in opts:
        if opt == '-h':
            usage_diffMaps()
            sys.exit(-1)
        elif opt in ("-i", "--inp1"):
            inp_file1 = arg
        elif opt in ("-j", "--inp2"):
            inp_file2 = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-t", "--type"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file1 = arg
        elif opt in ("-q", "--pdb2"):
            pdb_file2 = arg
        else:
            assert False, usage_diffMaps()

    # Input data matrix and PDB file are mandatory!
    if inp_file1 is None or inp_file2 is None or pdb_file1 is None:
        usage_diffMaps()
        sys.exit(-1)

    # Assign a default name if the user forgets the output file prefix.
    if out_file is None:
        out_file = "diff-map"

    # The user may prefer not to submit a title for the output.
    if sel_type is None:
        sel_type = "ndcc"

    return inp_file1, inp_file2, out_file, sel_type, pdb_file1, pdb_file2


def diffMapApp():
    """
    This app helps to plot difference maps.
    It doesnt' assume that the the matrix sizes are equal in both maps.
    Moreover, if you provide a second pdb file, it can match residue
    numbers and names of Calpha atoms in the pdb files.
    Please, beware that the program does not do any sequence alignment.
    Therefore, if one of the proteins contains a mutation, the app won't work.
    """
    inp_file1, inp_file2, out_file, sel_type, pdb_file1, pdb_file2 = handle_arguments_diffMaps()
    print(f"""
@> Running 'Difference Map App'

@> Input file   : {inp_file1}
@> Input file   : {inp_file2}
@> PDB file 1   : {pdb_file1}
@> PDB file 2   : {pdb_file2}
@> Data type    : {sel_type}
@> Output       : {out_file}""")

    selectedAtomSet1 = parsePDB(pdb_file1, subset='ca')

    if sel_type == 'nlmi':
        ccMatrix1 = convertLMIdata2Matrix(inp_file1, writeAllOutput=False)
        ccMatrix2 = convertLMIdata2Matrix(inp_file2, writeAllOutput=False)

        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif sel_type == 'absndcc':
        ccMatrix1 = np.absolute(np.loadtxt(inp_file1, dtype=float))
        ccMatrix2 = np.absolute(np.loadtxt(inp_file2, dtype=float))
        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif sel_type == 'ndcc':
        ccMatrix1 = np.loadtxt(inp_file1, dtype=float)
        ccMatrix2 = np.loadtxt(inp_file2, dtype=float)
        minColorBarLimit = -2
        maxColorBarLimit = 2

    else:
        print("Error: Unknown matrix type!")
        print("What is the type of your correlation map: nlmi, ndcc or absndcc?")
        sys.exit(-1)

    # One has to check if the lengths of two matrices match each other.
    # Otherwise, there is a huge problem here. We have to match corresponding
    # residues and then do the subtraction.
    if (len(ccMatrix1) == len(ccMatrix2)) and (len(ccMatrix1[0]) == len(ccMatrix2[0])):
        ccDiffMatrix = np.subtract(ccMatrix1, ccMatrix2)
        overallUniformDifferenceMap(ccMatrix1, ccMatrix2,
                                    minColorBarLimit, maxColorBarLimit,
                                    out_file, " ", selectedAtomSet1)
        # overallCorrelationMap(ccDiffMatrix, minColorBarLimit, maxColorBarLimit, \
        #                                                out_file, " ", selectedAtomSet1)
        ##########################################################################
        # Check number of chains. If there are multiple chains, plot inter and
        # intra chain correlations
        chains = Counter(selectedAtomSet1.getChids()).keys()
        saveMatrix = False
        plotChains = True
        if len(chains) > 1 and plotChains is True:
            intraChainCorrelationMaps(ccDiffMatrix,
                                      minColorBarLimit, maxColorBarLimit,
                                      out_file, " ",
                                      selectedAtomSet1, saveMatrix)
            interChainCorrelationMaps(ccDiffMatrix,
                                      minColorBarLimit, maxColorBarLimit,
                                      out_file, " ",
                                      selectedAtomSet1, saveMatrix)
    else:
        print("@> Warning: Sizes of two matrices are not equal!")

        if pdb_file2 is None:
            print("@> Warning: You have to specify at least two pdb files when")
            print("@>          matrix sizes are not identical!")
        else:
            selectedAtomSet2 = parsePDB(pdb_file2, subset='ca')
            overallNonUniformDifferenceMap(ccMatrix1, ccMatrix2,
                                           minColorBarLimit, maxColorBarLimit,
                                           out_file, " ",
                                           selectedAtomSet1, selectedAtomSet2)


if __name__ == '__main__':
    diffMapApp()

"""
This module help to plot difference maps.
It assumes the matrix sizes in both maps are equal
and number of Calpha atoms in the pdb file is in agreement with this number!
In fact, these assumptions have to be checked. 
As many scientists do, I am in a hurry and I'll just make it work until it 
breaks. 
"""
import os
import sys
import getopt

from collections import Counter, OrderedDict
import correlationPlus as cP
import numpy as np
from prody import parsePDB
import matplotlib.pyplot as plt


def usage_diffMaps():
    """
    Show how to use this program!
    """
    print("\nExample minimal usage:\n")
    print("python diffMap.py -i 4z90-cross-correlations.txt -j 4z91-cross-correlations.txt -p 4z90.pdb \n")
    print("Arguments: -i: The first file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)")
    print("           -j: The second file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)")
    print("           -p: PDB file of the protein. (Mandatory)")
    print("           -q: A second PDB file for the other conformation if residues numbers are not same in two conformations. (Optional)")
    print("           -t: It can be ndcc, lmi or absndcc (absolute values of ndcc). Default value is ndcc (Optional)")
    print("           -o: This will be your output file. Output figures are in png format. (Optional)\n\n")

def handle_arguments():
    opts, args = getopt.getopt(sys.argv[1:],"hi:j:o:t:p:q:",["inp1=", "inp2=", "out=", "type=", "pdb=", "pdb2="])
    inp_file1 = None
    inp_file2 = None
    pdb_file1 = None
    pdb_file2 = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:j:o:t:p:q:",["inp1=", "inp2=", "out=", "type=", "pdb=", "pdb2="])
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

    #Input data matrix and PDB file are mandatory!
    if inp_file1==None or inp_file2==None or pdb_file1==None:
        usage_diffMaps()
        sys.exit(-1)

    #Assign a default name if the user forgets the output file prefix.
    if (out_file == None):
        out_file = "diff-map"

    #The user may prefer not to submit a title for the output.
    if (sel_type == None):
        sel_type = "ndcc"

    return (inp_file1, inp_file2, out_file, sel_type, pdb_file1, pdb_file2)

def main():
    (inp_file1, inp_file2, out_file, sel_type, pdb_file1, pdb_file2) = handle_arguments()

    selectedAtomSet1 = parsePDB(pdb_file1, subset='ca')


    if (sel_type == 'lmi'):
        ccMatrix1 = cP.convertLMIdata2Matrix(inp_file1, writeAllOutput=False)
        ccMatrix2 = cP.convertLMIdata2Matrix(inp_file2, writeAllOutput=False)
        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif ((sel_type == 'absndcc')):
        ccMatrix1 = np.absolute(np.loadtxt(inp_file1, dtype=float))
        ccMatrix2 = np.absolute(np.loadtxt(inp_file2, dtype=float))
        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif ((sel_type == 'ndcc')):
        ccMatrix1 = np.loadtxt(inp_file1, dtype=float)
        ccMatrix2 = np.loadtxt(inp_file2, dtype=float)
        minColorBarLimit = -2
        maxColorBarLimit = 2

    else:
        print("Error: Unknown matrix type!")
        print("What is the type of your correlation map: lmi, ndcc or absndcc?")
        sys.exit(-1)

    #One has to check if the lengths of two matrices match each other. 
    #Otherwise, there is a huge problem here. We have to match corresponding
    #residues and then do the subtraction. 
    if((len(ccMatrix1) == len(ccMatrix2)) and \
        (len(ccMatrix1[0]) == len(ccMatrix2[0]))):
        ccDiffMatrix = np.subtract(ccMatrix1, ccMatrix2)
        cP.overallUniformDifferenceMap(ccMatrix1, ccMatrix2, minColorBarLimit, maxColorBarLimit, \
                                                                out_file, " ", selectedAtomSet1)
        #cP.overallCorrelationMap(ccDiffMatrix, minColorBarLimit, maxColorBarLimit, \
        #                                                out_file, " ", selectedAtomSet1)
        ##########################################################################
        #Check number of chains. If there are multiple chains, plot inter and 
        #intra chain correlations
        chains = Counter(selectedAtomSet1.getChids()).keys()
        saveMatrix = False
        plotChains = True
        if((len(chains)>1) & (plotChains == True)):
            cP.intraChainCorrelationMaps(ccDiffMatrix, minColorBarLimit, maxColorBarLimit,\
                                            out_file, " ", selectedAtomSet1, saveMatrix)
            cP.interChainCorrelationMaps(ccDiffMatrix, minColorBarLimit, maxColorBarLimit,\
                                            out_file, " ", selectedAtomSet1, saveMatrix )
    else:
        print("@> Warning: Sizes of two matrices are not equal!")
        
        if(pdb_file2 == None):
            print("@> Warning: You have to specify at least two pdb files when")
            print("@>          matrix sizes are not identical!")
        else:
            selectedAtomSet2 = parsePDB(pdb_file2, subset='ca')
            cP.overallNonUniformDifferenceMap(ccMatrix1, ccMatrix2, \
                                minColorBarLimit, maxColorBarLimit, \
                                out_file, " ", \
                                selectedAtomSet1, selectedAtomSet2)




if __name__ == "__main__":
    main()
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

import correlationPlus as cP
import numpy as np
from prody import parsePDB


def usage():
    """
    Show how to use this program!
    """
    print("\nExample minimal usage:\n")
    print("python diffMap.py -i 4z90-cross-correlations.txt -j 4z91-cross-correlations.txt -p 4z90.pdb \n")
    print("Arguments: -i: The first file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)")
    print("           -j: The second file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)")
    print("           -p: PDB file of the protein. (Mandatory)")
    print("           -t: It can be dcc, lmi or absdcc (absolute values of dcc). Default value is dcc (Optional)")
    print("           -o: This will be your output file. Output figures are in png format. (Optional)\n\n")

def handle_arguments():
    opts, args = getopt.getopt(sys.argv[1:],"hi:j:o:t:p:",["inp1=", "inp2=", "out=", "type=", "pdb="])
    inp_file1 = None
    inp_file2 = None
    pdb_file = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:j:o:t:p:",["inp1=", "inp1=", "out=", "type=", "pdb="])
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
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
            pdb_file = arg
        else:
            assert False, usage()

    #Input data matrix and PDB file are mandatory!
    if inp_file1==None or inp_file2==None or pdb_file==None:
        usage()
        sys.exit(-1)

    #Assign a default name if the user forgets the output file prefix.
    if (out_file == None):
        out_file = "diff-map"

    #The user may prefer not to submit a title for the output.
    if (sel_type == None):
        sel_type = "dcc"

    return (inp_file1, inp_file2, out_file, sel_type, pdb_file)

def main():
    (inp_file1, inp_file2, out_file, sel_type, pdb_file) = handle_arguments()

    selectedAtoms = parsePDB(pdb_file, subset='ca')


    if (sel_type == 'lmi'):
        ccMatrix1 = cP.convertLMIdata2Matrix(inp_file1, writeAllOutput=False)
        ccMatrix2 = cP.convertLMIdata2Matrix(inp_file2, writeAllOutput=False)
        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif ((sel_type == 'absdcc')):
        ccMatrix1 = np.absolute(np.loadtxt(inp_file1, dtype=float))
        ccMatrix2 = np.absolute(np.loadtxt(inp_file2, dtype=float))
        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif ((sel_type == 'dcc')):
        ccMatrix1 = np.loadtxt(inp_file1, dtype=float)
        ccMatrix2 = np.loadtxt(inp_file2, dtype=float)
        minColorBarLimit = -2
        maxColorBarLimit = 2

    else:
        print("Error: Unknown matrix type!")
        print("What is the type of your correlation map: lmi, dcc or absdcc?")
        sys.exit(-1)

    ccDiffMatrix = np.subtract(ccMatrix1, ccMatrix2)
    #cP.overallDifferenceMap(ccMatrix1, ccMatrix2, minColorBarLimit, maxColorBarLimit, \
    #                                                       out_file, " ", selectedAtoms)
    cP.overallCorrelationMap(ccDiffMatrix, minColorBarLimit, maxColorBarLimit, \
                                                    out_file, " ", selectedAtoms)
    ##########################################################################
    #Check number of chains. If there are multiple chains, plot inter and 
    #intra chain correlations
    chains = Counter(selectedAtoms.getChids()).keys()
    saveMatrix = False
    plotChains = True
    if((len(chains)>1) & (plotChains == True)):
        intraChainCorrelationMaps(ccDiffMatrix, minColorBarLimit, maxColorBarLimit,\
                                        out_file, " ", selectedAtoms, saveMatrix)
        interChainCorrelationMaps(ccDiffMatrix, minColorBarLimit, maxColorBarLimit,\
                                        out_file, " ", selectedAtoms, saveMatrix )




if __name__ == "__main__":
    main()
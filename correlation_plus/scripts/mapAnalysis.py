import sys
import getopt
from collections import Counter

import numpy as np
from prody import parsePDB

from correlation_plus.correlationPlus import overallCorrelationMap, convertLMIdata2Matrix, distanceDistribution
from correlation_plus.correlationPlus import intraChainCorrelationMaps, interChainCorrelationMaps
from correlation_plus.correlationPlus import filterCorrelationMapByDistance, projectCorrelationsOntoProteinVMD


def usage_mapAnalysisApp():
    """
    Show how to use this program!
    """
    print("""
Example usage:
mapAnalysis -i 4z90-cross-correlations.txt -p 4z90.pdb

Arguments: -i: A file containing normalized dynamical cross correlations in matrix format. (Mandatory)
           -p: PDB file of the protein. (Mandatory)
           -s: It can be ndcc, lmi or absndcc (absolute values of ndcc). Default value is ndcc (Optional)
           -o: This will be your output file. Output figures are in png format. (Optional)
""")


def handle_arguments_mapAnalysisApp():
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:s:p:", ["help", "inp=", "out=", "sel=", "pdb="])
    except getopt.GetoptError:
        usage_mapAnalysisApp()
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            usage_mapAnalysisApp()
            sys.exit(-1)
        elif opt in ("-i", "--inp"):
            inp_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-s", "--sel"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file = arg
        else:
            assert False, usage_mapAnalysisApp()

    # Input data matrix and PDB file are mandatory!
    if inp_file is None or pdb_file is None:
        usage_mapAnalysisApp()
        sys.exit(-1)

    # Assign a default name if the user forgets the output file prefix.
    if out_file is None:
        out_file = "correlation"

    # The user may prefer not to submit a title for the output.
    if sel_type is None:
        sel_type = "ndcc"

    return inp_file, out_file, sel_type, pdb_file


def mapAnalysisApp():
    print("@> Running 'Correlation Map App'")
    inp_file, out_file, sel_type, pdb_file = handle_arguments_mapAnalysisApp()
    print("\n@> Input file   :", inp_file)
    print("@> PDB file     :", pdb_file)
    print("@> Data type    :", sel_type)
    print("@> Output       :", out_file)

    ##########################################################################
    # Read PDB file
    # TODO: This is the only place where I use Prody.
    # Maybe, I can replace it with a library that only parses
    # PDB files. Prody does a lot more!
    selectedAtoms = parsePDB(pdb_file, subset='ca')

    ##########################################################################
    # Read data file and assign to a numpy array
    if sel_type == "dcc":
        ccMatrix = np.loadtxt(inp_file, dtype=float)
    elif sel_type == "absdcc":
        ccMatrix = np.absolute(np.loadtxt(inp_file, dtype=float))
    elif sel_type == "lmi":
        ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=True)
    else:
        print("Unknown matrix format!\n")
        sys.exit(-1)

    # Check the data type in the matrix.
    minCorrelationValue = np.min(ccMatrix)

    maxCorrelationValue = np.max(ccMatrix)

    if (minCorrelationValue < 0.0):
        # Assume that it is an nDCC file
        minColorBarLimit = -1.0
    else:
        # Assume that it is an LMI file
        minColorBarLimit = 0.0

    if (maxCorrelationValue > 1.0):
        print("This correlation map is not normalized!")
        # TODO: At this point, one can ask the user if s/he wants to normalize it!
        sys.exit(-1)
    ##########################################################################
    # Call overall correlation calculation
    maxColorBarLimit = 1.0
    overallCorrelationMap(ccMatrix, minColorBarLimit, maxColorBarLimit,
                          out_file, " ", selectedAtoms)

    plotDistributions = True
    if (plotDistributions):
        if (sel_type == "dcc"):
            distanceDistribution(ccMatrix, out_file, "nDCC", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)

        elif (sel_type == "absdcc"):
            distanceDistribution(ccMatrix, out_file, "Abs(nDCC)",
                                 selectedAtoms, absoluteValues=True, writeAllOutput=False)

        elif (sel_type == "lmi"):
            distanceDistribution(ccMatrix, out_file, "LMI", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)

        else:
            print("Warning: Unknows correlation data.\n")
            print("         Correlations can be ndcc, absndcc, lmi!\n")

    ##########################################################################
    # Check number of chains. If there are multiple chains, plot inter and
    # intra chain correlations
    chains = Counter(selectedAtoms.getChids()).keys()
    saveMatrix = False
    plotChains = True
    if (len(chains) > 1) and (plotChains is True):
        intraChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,
                                  out_file, " ", selectedAtoms, saveMatrix)
        interChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,
                                  out_file, " ", selectedAtoms, saveMatrix)

    # Here, we can filter some correlation values closer than a distance.
    # Typically, it is supposed to filter out the correlation within the
    # same secondary structure etc.
    filterByDistance = False
    if filterByDistance:
        distanceValue = 5.0
        ccMatrix = filterCorrelationMapByDistance(ccMatrix, out_file, " ",
                                                  selectedAtoms, distanceValue,
                                                  absoluteValues=False,
                                                  writeAllOutput=False)

    # Overall projection
    projectCorrelationsOntoProteinVMD(ccMatrix, out_file,
                                      selectedAtoms, valueFilter=0.75,
                                      absoluteValues=True, writeAllOutput=True)

    print("\n@> Program finished successfully!\n")


if __name__ == "__main__":
    mapAnalysisApp()

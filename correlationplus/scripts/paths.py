##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   dynamical correlations maps of proteins.                  #
# Authors: Mustafa Tekpinar                                                   #
# Copyright (C) Mustafa Tekpinar 2017-2018                                        #
# Copyright (C) CNRS-UMR3528, 2019                                                #
# Copyright (C) Institut Pasteur Paris, 2020-2021                                 #
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
from correlationplus.pathAnalysis import pathAnalysis, writePath2VMDFile
from correlationplus.pathAnalysis import mapResid2ResIndex
# from .visualizemap import handle_arguments_visualizemapApp


def usage_pathAnalysisApp():
    """
    Show how to use this program!
    """
    print("""
Example usage:
correlationplus paths -i 4z90-cross-correlations.txt -p 4z90.pdb -b A78 -e A230

Arguments: -i: A file containing dynamical correlations in 
               matrix format. (Mandatory)

           -p: PDB file of the protein. (Mandatory)
           
           -t: Type of the matrix. It can be ndcc, lmi or absndcc 
               (absolute values of ndcc). Default value is ndcc (Optional)

           -o: This will be your output file. Output figures are in png format. 
               (Optional)

           -v: Value filter. The values lower than this value in the map will be 
               considered as zero. Default is 0.3. (Optional)

           -d: Distance filter. The residues with distances higher than this value 
               will be considered as zero. Default is 7.0 Angstrom. (Optional)
                              
           -b: ChainID and residue ID of the beginning (source)  residue (Ex: A41). (Mandatory)

           -e: ChainID and residue ID of the end (sink/target) residue (Ex: B41). (Mandatory)

           -n: Number of shortest paths to write to tcl or pml files. Default is 1. (Optional)
""")


def handle_arguments_pathAnalysisApp():
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None
    val_fltr = None
    dis_fltr = None
    src_res = None
    trgt_res = None
    num_path = None

    try:
        opts, args = getopt.getopt(sys.argv[2:], "hi:o:t:p:v:d:b:e:n:", \
            ["help", "inp=", "out=", "type=", "pdb=", "value=", "distance", "beginning=", "end=", "distance=", "npaths"])
    except getopt.GetoptError:
        usage_pathAnalysisApp()
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            usage_pathAnalysisApp()
            sys.exit(-1)
        elif opt in ("-i", "--inp"):
            inp_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-t", "--type"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file = arg
        elif opt in ("-v", "--value"):
            val_fltr = arg
        elif opt in ("-d", "--distance"):
            dis_fltr = arg
        elif opt in ("-b", "--beginning"):
            src_res = arg
        elif opt in ("-e", "--end"):
            trgt_res = arg
        elif opt in ("-n", "--npaths"):
            num_path = arg
        else:
            assert False, usage_pathAnalysisApp()

    # Input data matrix and PDB file are mandatory!
    if inp_file is None or pdb_file is None:
        print("PDB file and a correlation matrix are mandatory!")
        usage_pathAnalysisApp()
        sys.exit(-1)

    # Assign a default name if the user forgets the output file prefix.
    if out_file is None:
        out_file = "paths"

    # The user may prefer not to submit a title for the output.
    if sel_type is None:
        sel_type = "ndcc"

    if val_fltr is None:
        val_fltr = 0.3

    if dis_fltr is None:
        dis_fltr = 7.0
    
    if num_path is None:
        num_path = 1

    if src_res is None:
        print("You have to specify a source resid!")
        usage_pathAnalysisApp()
        sys.exit(-1)

    if trgt_res is None:
        print("You have to specify a target resid!")
        usage_pathAnalysisApp()
        sys.exit(-1)

    return inp_file, out_file, sel_type, pdb_file, \
            val_fltr, dis_fltr, src_res, trgt_res, num_path


def pathAnalysisApp():
    inp_file, out_file, sel_type, pdb_file,val_fltr, \
    dis_fltr, src_res, trgt_res, num_paths\
            = handle_arguments_pathAnalysisApp()

    print(f"""
@> Running 'paths' app

@> Input file     : {inp_file}
@> PDB file       : {pdb_file}
@> Data type      : {sel_type}
@> Output         : {out_file}
@> Value filter   : {val_fltr}
@> Distance filter: {dis_fltr}
@> Source residue : {src_res}
@> Target residue : {trgt_res}
@> Number of paths: {num_paths}""")

    ##########################################################################
    # Read PDB file
    # TODO: This is the only place where I use Prody.
    # Maybe, I can replace it with a library that only parses
    # PDB files. Prody does a lot more!
    selectedAtoms = parsePDB(pdb_file, subset='ca')

    ##########################################################################
    # Read data file and assign to a numpy array
    if sel_type == "ndcc":
        ccMatrix = np.loadtxt(inp_file, dtype=float)
    elif sel_type == "absndcc":
        ccMatrix = np.absolute(np.loadtxt(inp_file, dtype=float))
    elif sel_type == "lmi":
        ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=True)
    else:
        print("Unknown data type: Type can only be ndcc, absndcc or lmi!\n")
        sys.exit(-1)

    sourceResid = src_res
    targetResid = trgt_res
    distanceMatrix = buildDistMatrix(selectedAtoms)
    resDict = mapResid2ResIndex(selectedAtoms)
    suboptimalPaths = pathAnalysis(ccMatrix, distanceMatrix, \
                                   val_fltr, dis_fltr,\
                                   resDict[sourceResid], resDict[targetResid], \
                                   selectedAtoms,\
                                   int(num_paths))

    out_file_full_name = out_file+"-source"+sourceResid+"-target"+targetResid+".vmd"
    writePath2VMDFile(suboptimalPaths, 
                    resDict[sourceResid], resDict[targetResid], \
                    pdb_file, out_file_full_name)
    
###############################################################################
# correlationPlus - Python module to plot dynamical correlations maps         #
#                   for proteins.                                             #
# Authors: Mustafa Tekpinar                                                   #
# Copyright Mustafa Tekpinar 2017-2018                                        #
# Copyright CNRS-UMR3528, 2019                                                #
# Copyright Institut Pasteur Paris, 2020                                       #
#                                                                             #
# This file is part of correlationPlus.                                       #
#                                                                             #
# correlationPlus is free software: you can redistribute it and/or modify     #
# it under the terms of the GNU Lesser General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# correlationPlus is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU LESSER General Public License for more details.                         #
#                                                                             #
# You should have received a copy of the GNU Lesser General Public License    #
# along with correlationPlus.  If not, see <https://www.gnu.org/licenses/>.   #
###############################################################################


import sys

import numpy as np
from prody import parsePDB

from correlationPlus.mapAnalysis import convertLMIdata2Matrix
from correlationPlus.centralityAnalysis import centralityAnalysis

from .mapAnalysis import handle_arguments_mapAnalysisApp


def centralityAnalysisApp():
    inp_file, out_file, sel_type, pdb_file = handle_arguments_mapAnalysisApp()

    print(f"""
@> Running 'CentralityAnalysis app

@> Input file   : {inp_file}
@> PDB file     : {pdb_file}
@> Data type    : {sel_type}
@> Output       : {out_file}""")

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

    valueFilter = 0.3
    centralityAnalysis(ccMatrix, valueFilter, out_file, "degree", selectedAtoms)
    centralityAnalysis(ccMatrix, valueFilter, out_file, "betweenness", selectedAtoms)
    centralityAnalysis(ccMatrix, valueFilter, out_file, "closeness", selectedAtoms)
    centralityAnalysis(ccMatrix, valueFilter, out_file, "current_flow_betweenness", selectedAtoms)
    centralityAnalysis(ccMatrix, valueFilter, out_file, "current_flow_closeness", selectedAtoms)
    centralityAnalysis(ccMatrix, valueFilter, out_file, "eigenvector", selectedAtoms)

##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   dynamical correlations maps of proteins.                  #
# Authors: Mustafa Tekpinar                                                   #
# Copyright Mustafa Tekpinar 2017-2018                                        #
# Copyright CNRS-UMR3528, 2019                                                #
# Copyright Institut Pasteur Paris, 2020                                      #
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

from prody import calcANM, calcGNM
from prody import calcCrossCorr
from numpy import savetxt

def calcENMnDCC(selectedAtoms, saveMatrix=True, out_file="nDCC.dat", method="ANM", nmodes=100):
    """
        Calculate normalized dynamical cross-correlations based on elastic 
        network model.
    
    Parameters
    ----------
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    saveMatrix: bool
        If True, an output file for the correlations will be written
        be written. 
    out_file: string
        Output file name for the data matrix. 
        Default value is nDCC.dat
    method: string
        This string can only take two values: ANM or GNM
        ANM us the default value.
    nmodes: int
        100 modes are default for normal mode based nDCC calculations.
    Returns
    -------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    """
    if(method == "ANM"):
        modes, sel = calcANM(selectedAtoms, n_modes=nmodes)
    elif (method == "GNM"):
        modes, sel = calcGNM(selectedAtoms, n_modes=nmodes)
 
    ccMatrix = calcCrossCorr(modes, n_cpu=1, norm=True)

    if(saveMatrix == True):
        savetxt(out_file, ccMatrix, fmt='%.6f')

    return ccMatrix




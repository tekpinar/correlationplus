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

import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numba

def calcENMnDCC(selectedAtoms, cut_off, method="ANM", nmodes=100, \
                normalized=True, saveMatrix=True, out_file="nDCC.dat"):
    """
        Calculate normalized dynamical cross-correlations based on elastic 
        network model.
    
    Parameters
    ----------
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    cut_off: int
        Cutoff radius in Angstrom unit for ANM or GNM. 
        Default value is 15 for ANM and 10 for GNM. 
    method: string
        This string can only take two values: ANM or GNM
        ANM us the default value.
    nmodes: int
        100 modes are default for normal mode based nDCC calculations.
    saveMatrix: bool
        If True, an output file for the correlations will be written
        be written. 
    out_file: string
        Output file name for the data matrix. 
        Default value is nDCC.dat
    Returns
    -------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    """
    if(method == "ANM"):
        modes, sel = calcANM(selectedAtoms, cutoff=cut_off, n_modes=nmodes)
    elif (method == "GNM"):
        modes, sel = calcGNM(selectedAtoms, cutoff=cut_off, n_modes=nmodes)
    
    if(normalized==True):
        ccMatrix = calcCrossCorr(modes, n_cpu=1, norm=True)
        out_file = "n"+out_file
    else:
        ccMatrix = calcCrossCorr(modes, n_cpu=1, norm=False)
    if(saveMatrix == True):
        np.savetxt(out_file, ccMatrix, fmt='%.6f')

    return ccMatrix

def calcMDnDCC(topology, trajectory, startingFrame=0, endingFrame=(-1),\
                  normalized=True, saveMatrix=True, out_file="DCC"):
    """
        Calculate normalized dynamical cross-correlations when a topology
        and a trajectory file is given. 
    Parameters
    ----------
    topology: string
        A PDB file.
    trajectory: string
        A trajectory file in dcd, xtc or trr format.
    startingFrame: int
        You can specify this value if you want to exclude some initial frames
        from your cross-correlation calculations. Default value is 0.
    endingFrame: int
        You can specify this value if you want to calculate cross-correlation 
        calculations up to a certain ending frame. Default value is -1 and it 
        indicates the last frame in your trajectory.
    normalized: bool
        Default value is True and it means that the cross-correlation matrix
        will be normalized. 
    saveMatrix: bool
        If True, cross-correlation matrix will be written to an output file. 
    out_file: string
        Output file name for the cross-correlation matrix. 
        Default value is DCC and the file extension is .dat. 
    Returns
    -------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    """
    #Create the universe (That sounds really fancy :)
    universe = mda.Universe(topology, trajectory)

    # Create an atomgroup from the alpha carbon selection
    calphas = universe.select_atoms("protein and name CA")
    N = calphas.n_atoms
    print("@> Parsed "+str(N)+" Calpha atoms.")
    # Set your frame window for your trajectory that you want to analyze
    #startingFrame = 0
    if(endingFrame==(-1)):
        endingFrame = universe.trajectory.n_frames
    skip = 1 
    
    Rvector = []
    # Iterate through the universe trajectory
    for timestep in universe.trajectory[startingFrame:endingFrame:skip]:
        Rvector.append(calphas.positions.flatten())
    ##############################################

    #Perform Calpha alignment
    print("@> Aligning only Calpha atoms to the initial frame!")
    alignment = align.AlignTraj(universe, universe, select="protein and name CA", in_memory=True)
    alignment.run()

    #I reassign this bc in ccMatrix calculation, we may skip some part of the trajectory!
    N_Frames = len(Rvector)

    R_average = np.mean(Rvector, axis=0)
    print("@> Calculating cross-correlation matrix!")
    ccMatrix = DCCmatrixCalculation(N, np.array(Rvector), R_average)

    #Do the averaging
    ccMatrix=ccMatrix/float(N_Frames)


    if(normalized==True):
        cc_normalized=np.zeros((N, N), np.double)
        for i in range(0, N):
            for j in range (i, N):
                cc_normalized[i][j]=ccMatrix[i][j]/(( (ccMatrix[i][i])*(ccMatrix[j][j]) )**0.5)
                cc_normalized[j][i]=cc_normalized[i][j]
        
        if(saveMatrix==True):
            np.savetxt("n"+out_file, cc_normalized, fmt='%.6f')

        return cc_normalized
    else:
        for i in range(0, N):
            for j in range (i+1, N):
                ccMatrix[j][i] = ccMatrix[i][j]
        if(saveMatrix==True):
            np.savetxt(out_file, ccMatrix, fmt='%.6f')
        return ccMatrix

@numba.njit
def DCCmatrixCalculation(N, Rvector, R_average):
    """
        This function calculates upper triangle of dynamical cross-correlation
        matrix. 
    """
    ccMatrix = np.zeros((N, N), np.double)
    for k in range(0, len(Rvector)):
        if(k%100==0):
            print("Frame: "+str(k))
        deltaR=np.subtract(Rvector[k], R_average)
        for i in range(0, N):
            ind_3i=3*i
            for j in range (i, N):
                ind_3j=3*j
                ccMatrix[i][j]+= (deltaR[ind_3i]*deltaR[ind_3j] + \
                            deltaR[ind_3i+1]*deltaR[ind_3j+1] +\
                            deltaR[ind_3i+2]*deltaR[ind_3j+2])
    return ccMatrix

##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   correlations maps of proteins.                            #
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
from prody import ANM, calcANM, calcGNM
from prody import calcCrossCorr

import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numba

def writeSparseCorrData(out_file, cMatrix, \
                        selectedAtoms, \
                        Ctype: bool, 
                        symmetric: bool):
    """
        This function writes correlation data in sparse format.

        In a sparse matrix, only nonzero elements of the matrix are 
        given in 3 columns format:
        i j C_ij
        i and j are the indices of the matrix positions (or residue indices,
        not residue IDs given in PDB files).
        It returns nothing. 

    Parameters
    ----------
    out_file: string
        Correlation file to write.
    cMatrix: A numpy square matrix of floats
        Correlation matrix.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    Ctype: boolean
        If Ctype=True, location indices i and j indices start from 0.
        Otherwise, it is assumed to be starting from 1.
    symmetric: boolean
        If you select it True, it will write the upper (or lower) triangle.

    Returns
    -------
    Nothing. 

    """
    data_file = open(out_file, 'w')
    n = selectedAtoms.numAtoms()
 
    if(symmetric == True):
        if (Ctype == True):
            for i in range(0, n):
                for j in range(i, n):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format(i, j, cMatrix[i][j]))
        else:
            for i in range(0, n):
                for j in range(i, n):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format((i+1), (j+1), cMatrix[i][j]))
    else:
        if (Ctype == True):
            for i in range(0, n):
                for j in range(0, n):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format(i, j, cMatrix[i][j]))
        else:
            for i in range(0, n):
                for j in range(0, n):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format((i+1), (j+1), cMatrix[i][j]))

    print("@> Finished writing correlation matrix in sparse format!")

    data_file.close()

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
    saveSparse = False
    if method == "ANM":
        modes, sel = calcANM(selectedAtoms, cutoff=cut_off, n_modes=nmodes)
    elif method == "GNM":
        modes, sel = calcGNM(selectedAtoms, cutoff=cut_off, n_modes=nmodes)
    else:
        print("@> ERROR: Unknown method! Select ANM or GNM!")
        sys.exit(-1)
    if normalized:
        ccMatrix = calcCrossCorr(modes, n_cpu=1, norm=True)
        out_file = "n" + out_file
    else:
        ccMatrix = calcCrossCorr(modes, n_cpu=1, norm=False)
    if saveMatrix:
        if(saveSparse):
            writeSparseCorrData(out_file, ccMatrix, \
                                selectedAtoms, Ctype=True,\
                                symmetric=True)
        else:
            np.savetxt(out_file, ccMatrix, fmt='%.6f')

    return ccMatrix


def calcMDnDCC(topology, trajectory, startingFrame=0, endingFrame=(-1),
               normalized=True, alignTrajectory=True,
               saveMatrix=True, out_file="DCC"):
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
    alignTrajectory: bool
        Default value is True and it means that all frames in the trajectory 
        will be aligned to the initial frame.  
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
    # Create the universe (That sounds really fancy :)
    universe = mda.Universe(topology, trajectory)

    # Create an atomgroup from the alpha carbon selection
    calphas = universe.select_atoms("protein and name CA")
    N = calphas.n_atoms
    print(f"@> Parsed {N} Calpha atoms.")
    # Set your frame window for your trajectory that you want to analyze
    # startingFrame = 0
    if endingFrame == -1:
        endingFrame = universe.trajectory.n_frames
    skip = 1 
    
    Rvector = []
    # Iterate through the universe trajectory
    for timestep in universe.trajectory[startingFrame:endingFrame:skip]:
        Rvector.append(calphas.positions.flatten())
    ##############################################

    # Perform Calpha alignment
    if alignTrajectory:
        print("@> Aligning only Calpha atoms to the initial frame!")
        alignment = align.AlignTraj(universe, universe, select="protein and name CA", in_memory=True)
        alignment.run()

    # I reassign this bc in ccMatrix calculation, we may skip some part of the trajectory!
    N_Frames = len(Rvector)

    R_average = np.mean(Rvector, axis=0)
    print("@> Calculating cross-correlation matrix:")
    ccMatrix = DCCmatrixCalculation(N, np.array(Rvector), R_average)

    # Do the averaging
    ccMatrix = ccMatrix / float(N_Frames)

    if normalized:
        cc_normalized = np.zeros((N, N), np.double)
        for i in range(0, N):
            for j in range (i, N):
                cc_normalized[i][j] = ccMatrix[i][j] / ((ccMatrix[i][i] * ccMatrix[j][j]) ** 0.5)
                cc_normalized[j][i] = cc_normalized[i][j]
        
        if saveMatrix:
            np.savetxt("n" + out_file, cc_normalized, fmt='%.6f')

        return cc_normalized
    else:
        for i in range(0, N):
            for j in range(i + 1, N):
                ccMatrix[j][i] = ccMatrix[i][j]
        if saveMatrix:
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
        if k % 100 == 0:
            print("@> Frame: " + str(k))
        deltaR = np.subtract(Rvector[k], R_average)
        for i in range(0, N):
            ind_3i = 3 * i
            for j in range(i, N):
                ind_3j = 3 * j
                ccMatrix[i][j] += (deltaR[ind_3i]*deltaR[ind_3j] +
                                   deltaR[ind_3i + 1] * deltaR[ind_3j + 1] +
                                   deltaR[ind_3i + 2] * deltaR[ind_3j + 2])
    return ccMatrix


def calcMD_LMI(topology, trajectory, startingFrame=0, endingFrame=(-1),
               normalized=True, alignTrajectory=True,
               saveMatrix=True, out_file="LMI"):
    """
        Calculate linear mutual information when a topology
        and a trajectory file is provided. 

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
        Default value is True and it means that the linear mutual information
        matrix will be normalized. 
    alignTrajectory: bool
        Default value is True and it means that all frames in the trajectory 
        will be aligned to the initial frame. 
    saveMatrix: bool
        If True, linear mutual information matrix will be written to an output
        file. 
    out_file: string
        Output file name for the linear mutual information matrix. 
        Default value is LMI and the file extension is .dat. 

    Returns
    -------
    lmiMatrix: A numpy square matrix of floats
        Linear mutual information matrix.

    """
    # Create the universe (That sounds really fancy :)
    universe = mda.Universe(topology, trajectory)

    # Create an atomgroup from the alpha carbon selection
    calphas = universe.select_atoms("protein and name CA")
    N = calphas.n_atoms
    print(f"@> Parsed {N} Calpha atoms.")
    # Set your frame window for your trajectory that you want to analyze
    #startingFrame = 0
    if endingFrame == -1:
        endingFrame = universe.trajectory.n_frames
    skip = 1 
    
    Rvector = []
    # Iterate through the universe trajectory
    for timestep in universe.trajectory[startingFrame:endingFrame:skip]:
        Rvector.append(calphas.positions.flatten())
    ##############################################

    if alignTrajectory:
        # Perform Calpha alignment
        print("@> Aligning only Calpha atoms to the initial frame!")
        alignment = align.AlignTraj(universe, universe,
                                    select="protein and name CA", in_memory=True)
        alignment.run()

    # I reassign this bc in lmiMatrix calculation, we may skip some part of the trajectory!
    N_Frames = len(Rvector)

    R_average = np.mean(Rvector, axis=0)
    print("@> Calculating linear mutual information matrix from MD trajectory:")

    ################### LMI matrix calculation part!    
    fullCovarMatrix = np.zeros((3 * N, 3 * N), np.double)
    for k in range(0, len(Rvector)):
        if k % 100 == 0:
            print("@> Frame: " + str(k))
        deltaR = np.subtract(Rvector[k], R_average)
        fullCovarMatrix += np.outer(deltaR, deltaR)

    # Do the averaging
    fullCovarMatrix = fullCovarMatrix / float(N_Frames)

    lmiMatrix = np.zeros((N, N), np.double)

    # Just to make the diagonal element 1.0 when normalized!
    np.fill_diagonal(lmiMatrix, 2000.0) 
    
    #TODO: This part can be put under @numba.njit to accelarate it!
    for i in range(0, N):
        ind_3i = 3 * i
        ind_3iPlus1 = 3 * (i + 1)
        for j in range(i + 1, N):
            ind_3j = 3 * j
            ind_3jPlus1 = 3 * (j + 1)

            # Diagonal element i
            subMatrix_C_i = fullCovarMatrix[ind_3i:ind_3iPlus1,
                                            ind_3i:ind_3iPlus1]
            detC_i = (np.linalg.det(subMatrix_C_i))

            # Diagonal element j
            subMatrix_C_j = fullCovarMatrix[ind_3j:ind_3jPlus1,
                                            ind_3j:ind_3jPlus1]
            detC_j = (np.linalg.det(subMatrix_C_j))

            # Copy matrix elements
            subMatrix_C_ij = np.zeros((6, 6), np.double)

            subMatrix_C_ij[0:3, 0:3] = subMatrix_C_i
            subMatrix_C_ij[3:6, 3:6] = subMatrix_C_j

            subMatrix_C_ij[0:3, 3:6] = fullCovarMatrix[ind_3i:ind_3iPlus1,
                                                       ind_3j:ind_3jPlus1]
            subMatrix_C_ij[3:6, 0:3] = np.transpose(subMatrix_C_ij[0:3, 3:6])

            detC_ij = (np.linalg.det(subMatrix_C_ij))
            
            lmiMatrix[i][j] = 0.5*(np.log(detC_i * detC_j / detC_ij))

            lmiMatrix[j][i] = lmiMatrix[i][j]
    #########################################################################

    if normalized:
        lmi_normalized = np.zeros((N, N), np.double)
        lmi_normalized = np.sqrt(1.0 - np.exp(-2.0 / 3.0 * lmiMatrix))
        
        if saveMatrix:
            np.savetxt("n" + out_file, lmi_normalized, fmt='%.6f')

        return lmi_normalized
    else:
        for i in range(0, N):
            for j in range(i + 1, N):
                lmiMatrix[j][i] = lmiMatrix[i][j]
        if saveMatrix:
            np.savetxt(out_file, lmiMatrix, fmt='%.6f')
        return lmiMatrix


def calcENM_LMI(selectedAtoms, cut_off, method="ANM", nmodes=100,
                normalized=True, saveMatrix=True, out_file="nDCC.dat"):
    """
        Calculate normalized linear mutual information matrix based on elastic 
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
    if method == "ANM":
        #modes, sel = calcANM(selectedAtoms, cutoff=cut_off, n_modes=nmodes)
        anm = ANM('ANM analysis')
        anm.buildHessian(selectedAtoms, cutoff=cut_off)
        anm.calcModes(n_modes=nmodes)

    # elif (method == "GNM"):
    #     modes, sel = calcGNM(selectedAtoms, cutoff=cut_off, n_modes=nmodes)
    else:
        print("@> ERROR: LMI can only be calculated from ANM modes!")
        sys.exit(-1)
      
    Rvector = []

    #######################################################################################
    N = selectedAtoms.numAtoms()
    scalingCoeff = 1.0

    for i in range(0, nmodes):
        newMode = (scalingCoeff / np.sqrt(anm[i].getEigval())) * anm[i].getEigvec()
        Rvector.append(selectedAtoms.getCoords().flatten() + newMode)
        #Rvector.append(selectedAtoms.getCoords()+np.reshape(newMode, (N, 3)))

    #######################################################################################
    # I reassign this bc in lmiMatrix calculation, we may skip some part of the trajectory!
    N_Frames = len(Rvector)

    R_average = np.mean(Rvector, axis=0)
    print("@> Calculating linear mutual information matrix from normal modes:")

    ################### LMI matrix calculation part!    
    fullCovarMatrix = np.zeros((3*N, 3*N), np.double)
    for k in range(0, len(Rvector)):
        if k % 99 == 0:
            print("@> Mode: " + str(k + 1))
        deltaR = np.subtract(Rvector[k], R_average)
        fullCovarMatrix += np.outer(deltaR, deltaR)

    # Do the averaging
    fullCovarMatrix = fullCovarMatrix / float(N_Frames)

    lmiMatrix = np.zeros((N, N), np.double)

    # Just to make the diagonal element 1.0 when normalized!
    np.fill_diagonal(lmiMatrix, 2000.0) 
    
    #TODO: This part can be put under @numba.njit to accelarate it!
    for i in range(0, N):
        ind_3i = 3 * i
        ind_3iPlus1 = 3 * (i + 1)
        for j in range(i + 1, N):
            ind_3j = 3 * j
            ind_3jPlus1 = 3 * (j + 1)

            # Diagonal element i
            subMatrix_C_i = fullCovarMatrix[ind_3i:ind_3iPlus1,
                                            ind_3i:ind_3iPlus1]
            detC_i = (np.linalg.det(subMatrix_C_i))

            # Diagonal element j
            subMatrix_C_j = fullCovarMatrix[ind_3j:ind_3jPlus1,
                                            ind_3j:ind_3jPlus1]
            detC_j = (np.linalg.det(subMatrix_C_j))

            # Copy matrix elements
            subMatrix_C_ij = np.zeros((6, 6), np.double)

            subMatrix_C_ij[0:3, 0:3] = subMatrix_C_i
            subMatrix_C_ij[3:6, 3:6] = subMatrix_C_j

            subMatrix_C_ij[0:3, 3:6] = fullCovarMatrix[ind_3i:ind_3iPlus1,
                                                       ind_3j:ind_3jPlus1]
            subMatrix_C_ij[3:6, 0:3] = np.transpose(subMatrix_C_ij[0:3, 3:6])

            detC_ij = (np.linalg.det(subMatrix_C_ij))
            
            lmiMatrix[i][j] = 0.5 * (np.log(detC_i * detC_j / detC_ij))

            lmiMatrix[j][i] = lmiMatrix[i][j]
    #########################################################################

    if normalized:
        lmi_normalized = np.zeros((N, N), np.double)
        lmi_normalized = np.sqrt(1.0 - np.exp(-2.0 / 3.0 * lmiMatrix))
        
        if saveMatrix:
            np.savetxt("n" + out_file, lmi_normalized, fmt='%.6f')

        return lmi_normalized
    else:
        for i in range(0, N):
            for j in range(i + 1, N):
                lmiMatrix[j][i] = lmiMatrix[i][j]
        if saveMatrix:
            np.savetxt(out_file, lmiMatrix, fmt='%.6f')
        return lmiMatrix

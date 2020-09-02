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




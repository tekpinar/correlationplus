##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   dynamical correlations maps of proteins.                  #
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

from prody import parsePDB

from correlationplus.calculate import *


def usage_calculateApp():
    """                                                                                                                                                                                       
    Show how to use this program!                                                                                                                                                             
    """
    print("""
Example usage: 

If you would like to calculate ANM-based normalized cross-correlations from a 
pdb file:                                                                                                                                                                                
correlationplus calculate -p 4z90.pdb

If you would like to calculate cross-correlations from a reference pdb file 
and a trajectory file (in dcd, xtc or trr formats):                                                                                                                                                                                
correlationplus calculate -p 4z90.pdb -f 4z90.xtc

If you would like to calculate ANM-based linear mutual information maps from a 
pdb file:                                                                                                                                                                                
correlationplus calculate -p 4z90.pdb -t lmi

If you would like to calculate a linear mutual information map from a reference pdb file 
and a trajectory file (in dcd, xtc or trr formats):                                                                                                                                                                                
correlationplus calculate -p 4z90.pdb -f 4z90.xtc -t lmi

Arguments:                                                                                         
           -p: PDB file of the protein. (Mandatory)
           -f: A trajectory file in dcd, xtc or trr format. (Optional)
           -b: Beginning frame in the trajectory to calculate the 
               correlation map (Valid only if you are using a trajectory).
           -e: Ending frame in the trajectory to calculate the 
               correlation map (Valid only if you are using a trajectory).
           -m: Elastic network model to calculate the correlations. 
               It can be ANM or GNM. Default is ANM. (Optional)
               (Valid only when you don't have a trajectory file.)
           -n: Number of non-zero modes, when ANM or GNM is used. 
               Default is 100. (Optional)
               (It can not exceed 3N-6 in ANM and N-1 in GNM, where N 
               is number of Calpha atoms.)
           -c: Cutoff radius in Angstrom for ANM or GNM. (Optional) 
               Default is 15 for ANM and 10 for GNM. 
           -t: Type of the correlation matrix. It can be ndcc or lmi.
               Default value is ndcc. (Optional)
           -o: This will be your output data file.
               Default is DCC.dat. (Optional)
""")


def handle_arguments_calculateApp():
    method = None
    out_file = None
    sel_type = None
    pdb_file = None
    trj_file = None
    beg_frm = 0
    end_frm = -1
    num_mod = 100
    cut_off = None

    try:
        opts, args = getopt.getopt(sys.argv[2:], "hm:o:t:p:f:b:e:n:c:",
                        ["help", "method=", "out=", "type=", "pdb=", "frames=", "beg=", "end=", "n_modes=", "cutoff="])
    except getopt.GetoptError:
        usage_calculateApp()
        print("@> ERROR: Unknown option encountered!")
        sys.exit(-1)

    for opt, arg in opts:
        if opt in ('-h', "--help"):
            usage_calculateApp()
            sys.exit(-1)
        elif opt in ("-m", "--method"):
            method = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-t", "--type"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file = arg
        elif opt in ("-f", "--frames"):
            trj_file = arg
        elif opt in ("-b", "--beg"):
            beg_frm = int(arg)
        elif opt in ("-e", "--end"):
            end_frm = int(arg)
        elif opt in ("-n", "--n_modes"):
            num_mod = int(arg)
        elif opt in ("-c", "--cutoff"):
            cut_off = int(arg)
        else:
            assert False, usage_calculateApp()

    # Input data matrix and PDB file are mandatory!                                                                                                                                           
    if pdb_file is None:
        print("\n@> ERROR: PDB file is mandatory!\n")
        usage_calculateApp()
        sys.exit(-1)

    # If the user has not provided a trajectory file:
    if trj_file is None:
        # It means that the user wants only ENM-based calculations.                                                                                                                       
        if method is None:
            method = "ANM"
        if (method != "ANM") and (method != "GNM"):
            print("@> ERROR: Unknown elastic network model!")
            print("@> ENM method can only be ANM or GNM!")
            sys.exit(-1)
    else:
        method = "MD"
        # Check the file extension.
        if not trj_file.lower().endswith(('.dcd', '.xtc', '.trr')):
            print("@> ERROR: Unrecognized trajectory type!")
            print("@> A trajectory file can only be in dcd, xtc or trr format!")
            usage_calculateApp()
            sys.exit(-1)

    # Assign a default matrix type                                                                                                                               
    if sel_type is None:
        sel_type = "ndcc"

    # Assign a default name if the user forgets the output file prefix.                                                                                                                      
    if out_file is None:
        if sel_type.lower() == "ndcc":
            out_file = "DCC"
        elif sel_type.lower() == "lmi":
            out_file = "LMI"
        else:
            print("@> ERROR: Unknown correlation matrix calculation requested!")
            print("@> This app can only calculate lmi or ndcc matrices.")
            print("@> Please check what you specified with -t option!")
            usage_calculateApp()
            sys.exit(-1)

    return method, out_file, sel_type, pdb_file, trj_file, beg_frm, end_frm, num_mod, cut_off


def calculateApp():
    method, out_file, sel_type, pdb_file, trj_file, beg_frm, end_frm, num_mod, cut_off = handle_arguments_calculateApp()
    print(f"""                                                                                                                                                                                
@> Running 'calculate App'
                                            
@> Method          : {method}
@> Output          : {out_file}
@> Correlation     : {sel_type}
@> PDB file        : {pdb_file}""")

    if out_file.lower().endswith(('.dat', '.txt')) is False:
        out_file = out_file + ".dat"

    if trj_file is None:
        print("@> Number of modes : " + str(num_mod))
        if cut_off is None and method == "ANM":
            cut_off = 15
        if cut_off is None and method == "GNM":
            cut_off = 10
        print("@> Cutoff radius   : " + str(cut_off))
        # Read pdb file
        selectedAtoms = parsePDB(pdb_file, subset='ca')
        if sel_type == "lmi":
            calcENM_LMI(selectedAtoms, cut_off,
                        method=method, 
                        nmodes=num_mod,
                        normalized=True,
                        saveMatrix=True,
                        out_file=out_file)
        else:
            calcENMnDCC(selectedAtoms, cut_off,
                        method=method, 
                        nmodes=num_mod,
                        normalized=True,
                        saveMatrix=True,
                        out_file=out_file)
        
    else:
        print("@> Trajectory file : " + trj_file)
        print("@> Beginning frame : " + str(beg_frm))
        print("@> Ending frame    : " + str(end_frm))
        if sel_type == "lmi":
            calcMD_LMI(pdb_file, trj_file,
                       startingFrame=beg_frm,
                       endingFrame=end_frm,
                       normalized=True,
                       alignTrajectory=True,
                       saveMatrix=True,
                       out_file=out_file)
        else:
            calcMDnDCC(pdb_file, trj_file,
                       startingFrame=beg_frm,
                       endingFrame=end_frm,
                       normalized=True,
                       alignTrajectory=True,
                       saveMatrix=True,
                       out_file=out_file)

    print("@> Correlation calculation finished successfully!")


if __name__ == "__main__":
    calculateApp()

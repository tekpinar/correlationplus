###############################################################################                                                                                                               
# correlationPlus - Python package to plot dynamical correlations maps        #                                                                                                               
#                   for proteins.                                             #                                                                                                               
# Authors: Mustafa Tekpinar                                                   #                                                                                                               
# Copyright Mustafa Tekpinar 2017-2018                                        #                                                                                                               
# Copyright CNRS-UMR3528, 2019                                                #                                                                                                               
# Copyright Institut Pasteur Paris, 2020                                      #                                                                                                              
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
import getopt

from prody import parsePDB

from correlationPlus.calculate import calcENMnDCC

def usage_calculateApp():
    """                                                                                                                                                                                       
    Show how to use this program!                                                                                                                                                             
    """
    print("""                                                                                                                                                                                 
Example for basic usage:                                                                                                                                                                                
correlationPlus calculate -p 4z90.pdb                                                                                                                        
                                                                                                                                                                                              
Arguments:                                                                                         
           -p: PDB file of the protein. (Mandatory)
           -m: Method to calculate correlations. It can be ANM or GNM (Default is ANM) . (Optional)                                                                                                                                           
           -t: Type of the matrix. It can be ndcc, lmi or absndcc (absolute values of ndcc). Default value is ndcc (Optional)                                                                 
           -o: This will be your output data file. Default is correlationMap.dat. (Optional)                                                                                                    
""")

def handle_arguments_calculateApp():
    method   = None
    out_file = None
    sel_type = None
    pdb_file = None

    try:
        opts, args = getopt.getopt(sys.argv[2:], "hm:o:t:p:", ["help", "method=", "out=", "type=", "pdb="])
    except getopt.GetoptError:
        usage_calculateApp()

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
        else:
            assert False, usage_calculateApp()

    # Input data matrix and PDB file are mandatory!                                                                                                                                           
    if pdb_file is None:
        usage_calculateApp()
        sys.exit(-1)

    # Assign a default name if the user forgets the output file prefix.                                                                                                                       
    if method is None:
        method = "ANM"
    if (method != "ANM") and (method != "GNM"):
        print("ERROR: Unknown elastic network model!")
        print("       ENM method can only be ANM or GNM!")
        sys.exit(-1)

    # Assign a default name if the user forgets the output file prefix.                                                                                                                       
    if out_file is None:
        out_file = "correlationMap.dat"

    # Assign a default matrix type                                                                                                                               
    if sel_type is None:
        sel_type = "ndcc"

    return method, out_file, sel_type, pdb_file
def calculateApp():
    method, out_file, sel_type, pdb_file = handle_arguments_calculateApp()
    print(f"""                                                                                                                                                                                
@> Running 'calculate App'
                                            
@> Method       : {method}
@> Output       : {out_file}
@> Data type    : {sel_type}
@> PDB File     : {pdb_file}
    """)

    #Read pdb file
    selectedAtoms = parsePDB(pdb_file, subset='ca')
    calcENMnDCC(selectedAtoms, saveMatrix=True, out_file=out_file, method=method, nmodes=100)
    print("@> Correlation calculation finished succesfully!")

if __name__ == "__main__":
    calculateApp()

##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   dynamical correlations maps of proteins.                  #
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
from math import fabs, log
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from prody import writePDB
import networkx as nx

from itertools import islice
from correlationplus.centralityAnalysis import buildDynamicsNetwork

def k_shortest_paths(G, source, target, k, weight=None):
    return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

def writePath2PMLFile(suboptimalPaths, selectedAtoms, source, target, pdb, outfile):
    """
    Produces PML output files for visualizing suboptimal paths with PyMol.

    This function writes a pml file containing the (suboptimal) paths
    between a source and a target residue. 
    The output files can be visualized with PyMol program by issuing the
    following command, if PyMol is in your path:
    pymol paths-sourceA145-targetB145.pml

    Parameters
    ----------
    suboptimalPaths: list of lists
        Each path is a list containing the indices of residues on the path.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    source: int
        This is the source residue index. Conversion from chainIDResID to 
        index is performed internally by mapResid2ResIndex() function. 
    target: int
        This is the target residue index. Conversion from chainIDResID to 
        index is performed internally by mapResid2ResIndex() function. 
    pdb: string
        This is a the name of the pdb file you submitted.
    out_file: string
        Prefix of the output file. Defaults value is 
        paths-source{chainIDResID}-target{chainIDResID}.pml

    Returns
    -------
    Nothing
    """

    pmlfile=open(outfile, "w+")
    pmlfile.write("from pymol.cgo import *\n")
    pmlfile.write("from pymol import cmd\n")

    pmlfile.write(f"load {pdb} \n\n")
    pmlfile.write("set cartoon_transparency, 0.5\n")
    pmlfile.write("cartoon type = tube\n")
    


    vdw_representation_string = "show spheres, chain {0:s} and resi {1:d} and name ca\n"
    draw_string = ("CYLINDER,  {0:.3f}, {1:.3f}, {2:.3f},\
    {3:.3f}, {4:.3f}, {5:.3f}, {6:.3f},\
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, \n ")
    endpoints_string = "select endpoints, \
                (chain {0:s} and resi {1:d} and name ca) \
             or (chain {2:s} and resi {3:d} and name ca)\n"
    pmlfile.write(endpoints_string.format(\
                                    selectedAtoms.getChids()[source],\
                                    selectedAtoms.getResnums()[source],\
                                    selectedAtoms.getChids()[target],\
                                    selectedAtoms.getResnums()[target]))
    pmlfile.write("set sphere_color, blue, endpoints\n")
    pmlfile.write("set sphere_scale, 1.25, endpoints\n")
    pmlfile.write(vdw_representation_string.\
                format(selectedAtoms.getChids()[source],\
                       selectedAtoms.getResnums()[source]))
    pmlfile.write(vdw_representation_string.\
                format(selectedAtoms.getChids()[target],\
                       selectedAtoms.getResnums()[target]))

    pmlfile.write("set sphere_scale, 0.75\n")
    pmlfile.write("set sphere_color, blue\n")
    k = 0
    for path in suboptimalPaths:
        for atom in path:
            pmlfile.write(vdw_representation_string.format(selectedAtoms.getChids()[atom],
                                                            selectedAtoms.getResnums()[atom]))

        pmlfile.write("python\n")
        # Draw cylinders
        pmlfile.write(f"path_{k+1} = [ \n")
        for i in range (0, (len(path)-1)):  
            pmlfile.write(draw_string.format(\
                        selectedAtoms.getCoords()[path[i]][0],
                        selectedAtoms.getCoords()[path[i]][1],
                        selectedAtoms.getCoords()[path[i]][2],
                        selectedAtoms.getCoords()[path[i+1]][0],
                        selectedAtoms.getCoords()[path[i+1]][1],
                        selectedAtoms.getCoords()[path[i+1]][2],0.5))
        pmlfile.write("]\n")
        pmlfile.write(f"cmd.load_cgo(path_{k+1},'path_{k+1}')\n")
        pmlfile.write(f"cmd.set(\"cgo_line_width\",2.0,'path_{k+1}')\n")
        pmlfile.write("python end\n")
        k = k+1
    pmlfile.close()

def writePath2VMDFile(suboptimalPaths, selectedAtoms, source, target, pdb, outfile):
    """
    Produces VMD output files for visualizing suboptimal paths.

    This function writes a tcl file containing the (suboptimal) paths
    between a source and a target residue. CA atoms on each path is 
    colored in a different color. 
    The output files can be visualized with VMD (Visual Molecular
    dynamics) program as follows.
    i) Load your pdb file, whether via command line or graphical interface.
    ii) Go to Extensions -> Tk Console and then
    iii) source path-sourceA41-targetB41.tcl
    It can take some time to load the general script.

    Parameters
    ----------
    suboptimalPaths: list of lists
        Each path is a list containing the indices of residues on the path.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    source: int
        This is the source residue index. Conversion from chainIDResID to 
        index is performed internally by mapResid2ResIndex() function. 
    target: int
        This is the target residue index. Conversion from chainIDResID to 
        index is performed internally by mapResid2ResIndex() function. 
    pdb: string
        This is a the name of the pdb file you submitted.
    out_file: string
        Prefix of the output file. Defaults value is 
        paths-source{chainIDResID}-target{chainIDResID}.tcl

    Returns
    -------
    Nothing
    """

    vmdfile=open(outfile, "w+")
    vmdfile.write(f"mol new {pdb}\n")
    vmdfile.write("mol delrep 0 top\n\n")

    vmdfile.write("#Set general view\n")
    vmdfile.write("mol representation Tube\n")
    vmdfile.write("mol color chain\n")
    vmdfile.write("mol selection {all}\n")
    vmdfile.write("mol material Transparent\n")
    vmdfile.write("mol addrep top\n\n")



    vdw_representation_string = "mol representation VDW 1.0 20\n" + \
                                "mol color ColorID 0\n" + \
                                "mol selection \"chain {0:s} and resid {1:d} and name CA\"\n" + \
                                "mol material Glossy\n"+\
                                "mol addrep top\n"
    draw_string = "draw cylinder " \
                  " [lindex [[atomselect top \"chain {0:s} and resid {1:d} and name CA\"] get {{x y z}}] 0]" \
                  " [lindex [[atomselect top \"chain {2:s} and resid {3:d} and name CA\"] get {{x y z}}] 0]" \
                  " radius {4:.3f}\n"
    vmdfile.write("#Set source atoms\n")
    vmdfile.write(vdw_representation_string.\
                format(selectedAtoms.getChids()[source],\
                       selectedAtoms.getResnums()[source]))

    vmdfile.write("#Set target atoms\n")
    vmdfile.write(vdw_representation_string.\
                format(selectedAtoms.getChids()[target],\
                       selectedAtoms.getResnums()[target]))

    k = 0
    for path in suboptimalPaths:
        vmdfile.write("#Set atoms on the path\n")
        vmdfile.write("mol representation VDW 0.6 20\n")
        vmdfile.write("mol color colorID "+str(k+7)+"\n")
        vmdfile.write("mol selection {(residue ")
        for atom in path:
            vmdfile.write(str(atom)+" ")
        vmdfile.write("and name CA)}\n")
        vmdfile.write("mol material Glossy\n")
        vmdfile.write("mol addrep top\n")

        for i in range (0, (len(path)-1)):
            vmdfile.write(draw_string.format(\
                        selectedAtoms.getChids()[path[i]],
                        selectedAtoms.getResnums()[path[i]],
                        selectedAtoms.getChids()[path[i+1]],
                        selectedAtoms.getResnums()[path[i+1]],0.5))
        
        k = k+1

    vmdfile.close()

def mapResid2ResIndex(selectedAtoms):
    """
    Map resid to residue index.

    The user enters a chainID and resid such as A41. This means the 
    Calpha atom of resid 41 in chainA. This function converts this value
    to a zero based index. 

    Parameters
    ----------
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.

    Return
    ------
    resDict: dictionary
        Basically, each element of the dictionary is like {'A41': 40,...}

    """

    # Build a dictionary to convert resid to residue indices.
    resnumList = selectedAtoms.getResnums()
    resindexList = selectedAtoms.getResindices()
    chidList = selectedAtoms.getChids()

    # Merge chainID and resid lists
    resnumList = list(map(str, resnumList))
    chainAndResid = [i + j for i, j in zip(chidList, resnumList)]

    # Create a dictionary for mathching 
    # residue numbers with residue indices
    resDict = dict(zip(chainAndResid, resindexList))

    return resDict

def pathAnalysis(graph, \
                 valueFilter, distanceFilter,\
                 sourceResid, targetResid, \
                 selectedAtoms, num_paths):
    """
    This function calculates paths between a source and target residue.

    This function calculates a path between a source (beginning) and a 
    target (end) residue from the graph build out of your ccMatrix.
    Parameters
    ----------
    graph: object
        It is a Networkx Graph object.
    valueFilter: float
        The ccMatrix values lower than the valueFilter will be ignored.
    distanceFilter: float
        The distance values higher than the distanceFilter will be ignored
        and they will not be considered as edges in a network. 
        This kind of value pruning may work for low conformational change MD
        simulations or ENM based calculations. However, if there are large
        scale structural changes, it will be necessary to eliminate the edges 
        based on contacts and their preservation in during the entire simulation. 
    sourceResid: str
        Chain and residue ID of the source residue. For example, A41, B145 etc..
    targetResid: str
        Chain and residue ID of the target residue. For example, A145, B41 etc.. 
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    num_paths: int
        Number of shortest paths to write to tcl or pml files. Minimum is 1.

    Returns
    -------
    Nothing

    """
    # # Build the graph
    # graph = buildDynamicsNetwork(ccMatrix, distanceMatrix, \
    #                    valueFilter, distanceFilter,\
    #                    selectedAtoms)
   

    suboptimalPaths = k_shortest_paths(graph, source=sourceResid, target=targetResid, k=num_paths, weight='weight')
    # #shortestPath = nx.shortest_path(graph, source=residue, target=targetResid, weight='weight', method='dijkstra')
    # shortestPathScoreList = []
    # for path in suboptimalPaths:
    #     shortestPathScore = 0.0
    #     for j in range(len(path)-1):
    #     #    print(path[j], end="\n")
    #     #    print(shortestPath[j+1], end="\n")
    #         shortestPathScore = shortestPathScore + ccMatrix[path[j]][path[j+1]]
    #     shortestPathScoreList.append(shortestPathScore)

    # spNP = np.array(shortestPath)
    # print(spNP)
    # print(shortestPathScoreList)
    # print("The shortest path score is: ", end='')
    # print(np.array(shortestPathScoreList).sum().round(3))
    k = 0
    for path in suboptimalPaths:
        k = k + 1
        path_length = nx.path_weight(graph, path, weight="weight")
        print("Path "+str(k)+" length: "+str(path_length))

    return suboptimalPaths


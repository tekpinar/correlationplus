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
from math import fabs, log
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from prody import writePDB
import networkx as nx

from itertools import islice
# def projectCentralitiesOntoProteinVMD(centrality, centralityArray, \
#                                     out_file, selectedAtoms, scalingFactor):
#     """
#     Produces VMD output files for visualizing protein centralities.

#     This function writes a tcl file and a PDB file that can be viewed in
#     VMD. Bfactor field of the protein contains the centrality information.
#     The first N residues with the highest centrality are highlighed in VDW
#     representation.  that  that contains the centralities on
#     on Bfactor field of the pdb.
#     The output files can be visualized with VMD (Visual Molecular
#     dynamics) program as follows.
#     i) Load your pdb file, whether via command line or graphical interface.
#     ii) Go to Extensions -> Tk Console and then
#     iii) source vmd-output-general.tcl
#     It can take some time to load the general script.

#     Parameters
#     ----------
#     centrality: string
#         It can have 'degree', 'betweenness', 'closeness',
#         'current_flow_betweenness' or 'current_flow_closeness'.
#     centralityArray: A numpy data array ?
#         It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
#     out_file: string
#         Prefix of the output file. According to the centralty measure, it will be
#         extended.
#     selectedAtoms: object
#         This is a prody.parsePDB object of typically CA atoms of a protein.
#     ScalingFactor: float
#         Sometimes, the values of the centrality arrays are too small.
#         The scaling factor multiplies the array to make the values visible in
#         the Bfactor colums.

#     Returns
#     -------
#     Nothing

#     """

#     # Write output in VMD format
#     # Writing the output is very important for further analyses such as
#     # inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
#     #
#     percentage = 0.10
#     numKeyResidues = int(percentage * len(selectedAtoms))
#     VMD_FILE = open(out_file + '_' + centrality + '.tcl', 'w')

#     VMD_FILE.write("mol new " + out_file + "_" + centrality + ".pdb" + "\n")
#     VMD_FILE.write("mol modstyle 0 0 Tube 0.5 25\n")
#     VMD_FILE.write("mol modcolor 0 0 Beta\n")
#     VMD_FILE.write("mol modmaterial 0 0 Glossy\n")

#     vdw_representation_string = "mol representation VDW 0.750000 25.000000\n" + \
#                                 "mol material Glossy\n" + \
#                                 "mol color Beta\n" + \
#                                 "mol selection \"chain {0:s} and resid {1:d} and name CA\"\n" + \
#                                 "mol addrep 0\n"
#     sortedList = np.flip(np.argsort(centralityArray))
#     for i in range(0, numKeyResidues):
#         # print(centralityArray[sortedList[i]])
#         VMD_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[sortedList[i]],
#                                                         selectedAtoms.getResnums()[sortedList[i]]))

#     selectedAtoms.setBetas([scalingFactor * i for i in centralityArray])

#     writePDB(out_file + '_' + centrality + '.pdb', selectedAtoms)

#     VMD_FILE.close()

# def projectCentralitiesOntoProteinPyMol(centrality, centralityArray, out_file, \
#                                         selectedAtoms, scalingFactor):
#     """
#     Produces PyMol output files for visualizing protein centralities.

#     This function writes a pml file and a PDB file that can be viewed in
#     VMD. Bfactor field of the protein contains the centrality information.
#     The first N residues with the highest centrality are highlighed in VDW
#     representation.  that  that contains the centralities on
#     on Bfactor field of the pdb.
#     The output files can be visualized with VMD (Visual Molecular
#     dynamics) program as follows: pymol output.pml

#     Parameters
#     ----------
#     centrality: string
#         It can have 'degree', 'betweenness', 'closeness',
#         'current_flow_betweenness' or 'current_flow_closeness'.
#     centralityArray: A numpy data array ?
#         It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
#     out_file: string
#         Prefix of the output file. According to the centralty measure, it will be
#         extended.
#     selectedAtoms: object
#         This is a prody.parsePDB object of typically CA atoms of a protein.
#     ScalingFactor: float
#         Sometimes, the values of the centrality arrays are too small.
#         The scaling factor multiplies the array to make the values visible in
#         the Bfactor colums.

#     Returns
#     -------
#     Nothing

#     """

#     percentage = 0.10
#     numKeyResidues = int(percentage * len(selectedAtoms))
#     PML_FILE = open(out_file + '_' + centrality + '.pml', 'w')
    
#     PML_FILE.write("load " + out_file + "_" + centrality + ".pdb" + "\n")
#     PML_FILE.write("cartoon type = tube\n")
#     PML_FILE.write("spectrum b\n")
#     PML_FILE.write("set sphere_scale, 0.75\n\n")

#     vdw_representation_string = "show spheres, chain {0:s} and resi {1:d} and name ca\n"

#     sortedList = np.flip(np.argsort(centralityArray))
#     for i in range(0, numKeyResidues):
#         # print(centralityArray[sortedList[i]])
#         PML_FILE.write(vdw_representation_string.\
#             format(selectedAtoms.getChids()[sortedList[i]],
#                     selectedAtoms.getResnums()[sortedList[i]]))

#     selectedAtoms.setBetas([scalingFactor * i for i in centralityArray])

#     writePDB(out_file + '_' + centrality + '.pdb', selectedAtoms)

#     PML_FILE.close()


# def plotCentralities(centrality, centralityArray, out_file, selectedAtoms, \
#                     scalingFactor):
#     """
#     Plots the centrality values on a 2D graph.

#     The centrality values are plotted on a 2D png file.
#     If there are at least two chains, the function produces a figure
#     for each chain.

#     Parameters
#     ----------
#     centrality: string
#         It can have 'degree', 'betweenness', 'closeness',
#         'current_flow_betweenness' or 'current_flow_closeness'.
#     centralityArray: A numpy data array ?
#         It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
#     out_file: string
#         Prefix of the output file. According to the centralty measure, it will be
#         extended.
#     selectedAtoms: object
#         This is a prody.parsePDB object of typically CA atoms of a protein.
#     ScalingFactor: float
#         Sometimes, the values of the centrality arrays are too small.
#         The scaling factor multiplies the array to make the values visible in
#         the Bfactor colums.

#     Return
#     ------
#     Nothing

#     """

#     chains = Counter(selectedAtoms.getChids()).keys()

#     if len(chains) > 1:
#         # Intra-chain
#         for chain in chains:
#             x = []
#             y = []
#             dst_file = f"{out_file}_{centrality}_chain{chain}"
#             for i in range(0, len(selectedAtoms.getResnums())):
#                 if selectedAtoms.getChids()[i] == chain:
#                     x.append(selectedAtoms[i].getResnum())
#                     y.append(centralityArray[i])

#             fig, ax = plt.subplots()
#             plt.title('Chain ' + chain)
#             plt.locator_params(axis='y', nbins=5)

#             plt.xticks(fontsize=16)
#             plt.yticks(fontsize=16)
#             ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
#             plt.ylabel(centrality.replace('_', ' '), fontsize=20)
#             plt.xlabel("Residue Number", fontsize=20)

#             # plt.plot(x, centralityArray '.', color='k')
#             plt.bar(x, y, color='k')
#             plt.tight_layout()
#             # plt.show()
#             plt.savefig(dst_file + '.png')
#             plt.close('all')
#     else:
#         dst_file = f"{out_file}_{centrality}"
#         fig, ax = plt.subplots()
#         # plt.locator_params(axis='y', nbins=4)

#         plt.xticks(fontsize=16)
#         plt.yticks(fontsize=16)
#         ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
#         plt.ylabel(centrality.replace('_', ' '), fontsize=20)
#         plt.xlabel("Residue Number", fontsize=20)

#         x = selectedAtoms.getResnums()
#         # plt.plot(x, centralityArray '.', color='k')
#         plt.bar(x, centralityArray, color='k')
#         plt.tight_layout()
#         # plt.xlim(xmin=0) #xmin will be removed in matplotlib 3.2
#         # plt.xlim(left=0)
#         # plt.show()
#         plt.savefig(dst_file + '.png')
#         plt.close('all')


# def projectCommunitiesOntoProteinVMD(sortedCommunities, out_file, selectedAtoms):
#     """
#     Produces VMD output files for visualizing protein communities.

#     This function writes a tcl file and a PDB file that can be viewed in
#     VMD. Occupancy field of the protein contains the community information.

#     The output files can be visualized with VMD (Visual Molecular
#     dynamics) program as follows.
#     i) Load your pdb file, whether via command line or graphical interface.
#     ii) Go to Extensions -> Tk Console and then
#     iii) source vmd-output-general.tcl

#     Parameters
#     ----------
#     sortedCommunities: Iterator over tuples of sets of nodes
#         It is a tuple of lists. Each list contain a community.
#     out_file: string
#         Prefix of the output file. According to the centralty measure, it will be
#         extended.
#     selectedAtoms: object
#         This is a prody.parsePDB object of typically CA atoms of a protein.
#     Returns
#     -------
#     Nothing

#     """
#     VMD_FILE = open(out_file + '_communities.tcl', 'w')

#     VMD_FILE.write("mol new " + out_file + '_communities.pdb' + "\n")
#     VMD_FILE.write("mol modstyle 0 0 Tube 0.5 25\n")
#     VMD_FILE.write("mol modcolor 0 0 Occupancy\n")
#     VMD_FILE.write("mol modmaterial 0 0 Glossy\n")

#     selectedAtoms.setOccupancies(0)
#     i=0
#     for item in sortedCommunities:
#         #print(item)
#         for residx in item:
#             selectedAtoms[residx].setOccupancy(i)
#         i = i + 1  

#     writePDB(out_file + '_communities.pdb', selectedAtoms)

#     VMD_FILE.close()

# def projectCommunitiesOntoProteinPyMol(sortedCommunities, out_file, selectedAtoms):
#     """
#     Produces PyMol output files for visualizing protein communities.

#     This function writes a pml file and a PDB file that can be viewed in
#     PyMol. Occupancy field of the protein contains the community information.

#     The output files can be visualized with PyMol program as follows:
#     pymol outputfile.pml

#     Parameters
#     ----------
#     sortedCommunities: Iterator over tuples of sets of nodes
#         It is a tuple of lists. Each list contain a community.
#     out_file: string
#         Prefix of the output file. According to the centralty measure, it will be
#         extended.
#     selectedAtoms: object
#         This is a prody.parsePDB object of typically CA atoms of a protein.
#     Returns
#     -------
#     Nothing

#     """
#     PML_FILE = open(out_file + '_communities.pml', 'w')

#     PML_FILE.write("load " + out_file + '_communities.pdb' + "\n")
#     PML_FILE.write("cartoon type = tube\n")
#     PML_FILE.write("spectrum q\n")
#     PML_FILE.write("set sphere_scale, 0.75\n\n")

#     selectedAtoms.setOccupancies(0)
#     i=0
#     for item in sortedCommunities:
#         #print(item)
#         for residx in item:
#             selectedAtoms[residx].setOccupancy(i)
#         i = i + 1  

#     writePDB(out_file + '_communities.pdb', selectedAtoms)

#     PML_FILE.close()
def k_shortest_paths(G, source, target, k, weight=None):
    return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

def buildDynamicsGraph(ccMatrix, distanceMatrix, \
                       valueFilter, distanceFilter,\
                       selectedAtoms):
    """
    This function calculates various network (graph) centralities of a protein.

    This function calculates some network centrality measures such as
    degree, betweenness, closeness, current flow betweenness and eigenvector.
    This function needs Python 3.6 or later to maintain dictionary order.!!!

    Parameters
    ----------
    ccMatrix: Numpy matrix
        It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
    distanceMatrix: Numpy matrix
        The distances between Calpha atoms of the protein stored in a matrix.
    valueFilter: float
        The ccMatrix values lower than the valueFilter will be ignored.
    distanceFilter: float
        The distance values higher than the distanceFilter will be ignored
        and they will not be considered as edges in a network. 
        This kind of value pruning may work for low conformational change MD
        simulations or ENM based calculations. However, if there are large
        scale structural changes, it will be necessary to eliminate the edges 
        based on contacts and their preservation in during the entire simulation.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.

    Returns
    -------
    A networkx graph object

    """
    # Create your  graph
    dynNetwork = nx.Graph()

    n = selectedAtoms.numAtoms()

    # Add all CA atoms as nodes
    for i in range(n):
        dynNetwork.add_node(i)

    # Add all pairwise interactions greater than the valueFilter as edges.
    # In addition, add only edges which has a distance of lower than the 
    # distance filter
    for i in range(n):
        for j in range(n):
            if fabs(ccMatrix[i][j]) > valueFilter and distanceMatrix[i][j] <= distanceFilter:
                dynNetwork.add_edge(i, j, weight=-log(fabs(ccMatrix[i][j])))
                # dynNetwork.add_edge(i, j, weight=fabs(correlationArray[i][j]))

    return dynNetwork

def writePath2PMLFile(suboptimalPaths, selectedAtoms, source, target, pdb, outfile):
    """
    Produces PML output files for visualizing suboptimal paths.

    This function writes a pml file containing the (suboptimal) paths
    between a source and a target residue. CA atoms on each path is 
    colored in a different color. 
    The output files can be visualized with PyMol program.

    Parameters
    ----------
    suboptimalPaths: list of lists
        Each path is a list containing the indices of residues on the path.
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
        pmlfile.write(f"path_{k} = [ \n")
        for i in range (0, (len(path)-1)):  
            pmlfile.write(draw_string.format(\
                        selectedAtoms.getCoords()[path[i]][0],
                        selectedAtoms.getCoords()[path[i]][1],
                        selectedAtoms.getCoords()[path[i]][2],
                        selectedAtoms.getCoords()[path[i+1]][0],
                        selectedAtoms.getCoords()[path[i+1]][1],
                        selectedAtoms.getCoords()[path[i+1]][2],0.5))
        pmlfile.write("]\n")
        pmlfile.write(f"cmd.load_cgo(path_{k},'path_{k}')\n")
        pmlfile.write(f"cmd.set(\"cgo_line_width\",2.0,'path_{k}')\n")
        pmlfile.write("python end\n")
        k = k+1
    pmlfile.close()

def writePath2VMDFile(suboptimalPaths, source, target, pdb, outfile):
    """
    Produces VMD output files for visualizing suboptimal paths.

    This function writes a tcl file containing the (suboptimal) paths
    between a source and a target residue. CA atoms on each path is 
    colored in a different color. 
    The output files can be visualized with VMD (Visual Molecular
    dynamics) program as follows.
    i) Load your pdb file, whether via command line or graphical interface.
    ii) Go to Extensions -> Tk Console and then
    iii) source vmd-output-general.tcl
    It can take some time to load the general script.

    Parameters
    ----------
    suboptimalPaths: list of lists
        Each path is a list containing the indices of residues on the path.
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
        paths-source{chainIDResID}-target{chainIDResID}.vmd

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

    vmdfile.write("#Set source atoms\n")
    vmdfile.write("mol representation VDW 1.0 20\n")
    vmdfile.write("mol color colorID 0\n")
    vmdfile.write("mol selection {(residue "+str(source)+" and name CA)}\n")
    vmdfile.write("mol material Glossy\n")
    vmdfile.write("mol addrep top\n\n")

    vmdfile.write("#Set target atoms\n")
    vmdfile.write("mol representation VDW 1.0 20\n")
    vmdfile.write("mol color colorID 0\n")
    vmdfile.write("mol selection {(residue "+str(target)+" and name CA)}\n")
    vmdfile.write("mol material Glossy\n")
    vmdfile.write("mol addrep top\n\n")

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
            vmdfile.write("draw cylinder [lindex [[atomselect top \"residue "+\
                str(path[i])+" and name CA\"] get {x y z}] 0] [lindex [[atomselect top \"residue "+\
                str(path[i+1])+" and name CA\"] get {x y z}] 0] radius 0.5 resolution 10 filled 0\n")
        
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

def pathAnalysis(ccMatrix, distanceMatrix,\
                 valueFilter, distanceFilter,\
                 sourceResid, targetResid, \
                 selectedAtoms, num_paths):
    """
    This function calculates paths between a source and target residue.

    This function calculates a path between a source (beginning) and a 
    target (end) residue from the graph build out of your ccMatrix.
    Parameters
    ----------
    ccMatrix: Numpy matrix
        It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
    distanceMatrix: Numpy matrix
        The distances between Calpha atoms of the protein stored in a matrix.
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
    # Build the graph
    dynNetwork = buildDynamicsGraph(ccMatrix, distanceMatrix, \
                       valueFilter, distanceFilter,\
                       selectedAtoms)
   
    shortestPathScoreList = []
    suboptimalPaths = k_shortest_paths(dynNetwork, source=sourceResid, target=targetResid, k=num_paths, weight='weight')
    #shortestPath = nx.shortest_path(dynNetwork, source=residue, target=targetResid, weight='weight', method='dijkstra')
    for path in suboptimalPaths:
        shortestPathScore = 0.0
        for j in range(len(path)-1):
        #    print(path[j], end="\n")
        #    print(shortestPath[j+1], end="\n")
            shortestPathScore = shortestPathScore + ccMatrix[path[j]][path[j+1]]
        shortestPathScoreList.append(shortestPathScore)

    # spNP = np.array(shortestPath)
    # print(spNP)
    # print(shortestPathScoreList)
    # print("The shortest path score is: ", end='')
    # print(np.array(shortestPathScoreList).sum().round(3))
    for path in suboptimalPaths:
        path_length = nx.path_weight(dynNetwork, path, weight="weight")
        print("Path length: "+str(path_length))

    return suboptimalPaths


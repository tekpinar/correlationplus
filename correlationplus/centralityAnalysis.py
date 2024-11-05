##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   correlation maps of proteins.                             #
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


def projectCentralitiesOntoProteinVMD(centrality, centralityArray, \
                                    out_file, selectedAtoms, scalingFactor):
    """
    Produces VMD output files for visualizing protein centralities.

    This function writes a tcl file and a PDB file that can be viewed in
    VMD. Bfactor field of the protein contains the centrality information.
    The first N residues with the highest centrality are highlighed in VDW
    representation.  that  that contains the centralities on
    on Bfactor field of the pdb.
    The output files can be visualized with VMD (Visual Molecular
    dynamics) program as follows.
    i) Load your pdb file, whether via command line or graphical interface.
    ii) Go to Extensions -> Tk Console and then
    iii) source vmd-output-general.tcl
    It can take some time to load the general script.

    Parameters
    ----------
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness' or 'current_flow_closeness'.
    centralityArray: A Python list.
        It is a floating point Python list of the centrality values.
    out_file: string
        Prefix of the output file. According to the centrality measure, it will be
        extended.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    ScalingFactor: float
        Sometimes, the values of the centrality arrays are too small.
        The scaling factor multiplies the array to make the values visible in
        the Bfactor colums.

    Returns
    -------
    Nothing

    """

    # Write output in VMD format
    # Writing the output is very important for further analyses such as
    # inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    #
    percentage = 0.10
    numKeyResidues = int(percentage * len(selectedAtoms))
    VMD_FILE = open(out_file + '_' + centrality + '.tcl', 'w')

    VMD_FILE.write("mol new " + out_file + "_" + centrality + ".pdb" + "\n")
    VMD_FILE.write("mol modstyle 0 0 Tube 0.5 25\n")
    VMD_FILE.write("mol modcolor 0 0 Beta\n")
    VMD_FILE.write("mol modmaterial 0 0 Glossy\n")

    vdw_representation_string = "mol representation VDW 0.750000 25.000000\n" + \
                                "mol material Glossy\n" + \
                                "mol color Beta\n" + \
                                "mol selection \"chain {0:s} and resid {1:d} and name CA\"\n" + \
                                "mol addrep 0\n"
    sortedList = np.flip(np.argsort(centralityArray))
    for i in range(0, numKeyResidues):
        # print(centralityArray[sortedList[i]])
        VMD_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[sortedList[i]],
                                                        selectedAtoms.getResnums()[sortedList[i]]))

    selectedAtoms.setBetas([scalingFactor * i for i in centralityArray])

    writePDB(out_file + '_' + centrality + '.pdb', selectedAtoms)

    VMD_FILE.close()

def projectCentralitiesOntoProteinPyMol(centrality, centralityArray, out_file, \
                                        selectedAtoms, scalingFactor):
    """
    Produces PyMol output files for visualizing protein centralities.

    This function writes a pml file and a PDB file that can be viewed in
    PyMol. Bfactor field of the protein contains the centrality information.
    The first N residues with the highest centrality are highlighed in VDW
    representation. 
    The output files can be visualized with PyMol program as follows: 
        pymol output.pml

    Parameters
    ----------
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness', 'current_flow_closeness', or 'eigenvector'.
    centralityArray: A numpy data array ?
        It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    ScalingFactor: float
        Sometimes, the values of the centrality arrays are too small.
        The scaling factor multiplies the array to make the values visible in
        the Bfactor colums.

    Returns
    -------
    Nothing

    """

    percentage = 0.10
    numKeyResidues = int(percentage * len(selectedAtoms))
    PML_FILE = open(out_file + '_' + centrality + '.pml', 'w')
    
    PML_FILE.write("load " + out_file + "_" + centrality + ".pdb" + "\n")
    PML_FILE.write("cartoon type = tube\n")
    PML_FILE.write("spectrum b\n")
    PML_FILE.write("set sphere_scale, 0.75\n\n")

    vdw_representation_string = "show spheres, chain {0:s} and resi {1:d} and name ca\n"

    sortedList = np.flip(np.argsort(centralityArray))

    for i in range(0, numKeyResidues):
        # print(centralityArray[sortedList[i]])
        PML_FILE.write(vdw_representation_string.\
            format(selectedAtoms.getChids()[sortedList[i]],
                    selectedAtoms.getResnums()[sortedList[i]]))

    PML_FILE.write("bg_color white\n")

    PML_FILE.write("png "+out_file + "_" + centrality + "_projection.png, dpi=300, ray=1\n")
    PML_FILE.close()
    selectedAtoms.setBetas([scalingFactor * i for i in centralityArray])

    writePDB(out_file + '_' + centrality + '.pdb', selectedAtoms)



def autoScaleCentralities(centralityArray,\
                          selectedAtoms, global_or_local):
    """
    Scales the centralities automatically.

    The centralities are scaled to [0, 1] according to the 
    c_scaled = c_i - c_min / (c_max - c_min) formula.

    Please note that the centrality quantity must always have 
    positive values for this function. 

    Parameters
    ----------
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness', 'current_flow_closeness' or 'eigenvector'.
    centralityArray: A Python list.
        It is a floating point Python list of the centrality values.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    global_or_local: string
        This string can only have one of the two strings: 'global', 'local'.
        If it is 'global', the scaling will be performed for the entire protein.
        This is particularly good for monomeric structures. 
        If 'local' is selected, the scaling will be performed for each chain 
        separateyl. This can be particularly good for multimeric structures. 

    Returns
    -------
    centralityArrayUpdated: A numpy data array

    """
    
    chains = Counter(selectedAtoms.getChids()).keys()

    idx = 0
    if (len(chains) > 1):
        if(global_or_local == 'global'):
            c_max = np.amax(centralityArray)
            # print (c_max)
            c_min = np.amin(centralityArray)
            # print (c_min)
            diff = (c_max - c_min)
            for i in range (len(centralityArray)):
                centralityArray[i] = (centralityArray[i] - c_min) / diff
            # print(centralityArray)
        elif(global_or_local == 'local'):
            for ch in chains:
                # Find beginning and end indices of the chains. 
                localCentralityArray = np.where(selectedAtoms.getChids() == ch)
                # print(localCentralityArray[0])
                beg = localCentralityArray[0][0]
                end = localCentralityArray[0][-1]
                # print(beg, end)
                # print(centralityArray[beg:end+1])
                c_max = np.amax(centralityArray[beg:end+1])
                c_min = np.amin(centralityArray[beg:end+1])
                diff = (c_max - c_min)

                for j in range(beg, end+1):
                    centralityArray[j] = (centralityArray[j] - c_min) / diff
        else:
            print("\n@> ERROR: Unknown global_or_local value!\n")
            print("@>        It can only be 'global' or 'local'")
    elif len(chains) == 1:
        c_max = np.amax(centralityArray)
        # print (c_max)
        c_min = np.amin(centralityArray)
        # print (c_min)
        diff = (c_max - c_min)
        for i in range (len(centralityArray)):
            centralityArray[i] = (centralityArray[i] - c_min) / diff
        # print(centralityArray)
    else:
        print("@> ERROR: Shame on you! There are not chain IDs in the pdb!")
        sys.exit(-1)
    return centralityArray


def projectCentralitiesOntoProteinPyMol_v2(centrality, centralityArray, out_file, \
                                            selectedAtoms, percentage=0.10):
    """
    Produces PyMol output files for visualizing protein centralities
    projected onto protein B factors column.

    This function writes a pml file and a PDB file that can be viewed in
    PyMol. Bfactor field of the protein contains the centrality information.
    The first N residues with the highest centrality are highlighed in VDW
    representation. 
    The output files can be visualized with PyMol program as follows: 
        pymol output.pml

    Parameters
    ----------
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness', 'current_flow_closeness' or 'eigenvector'.
    centralityArray: A Python list.
        It is a floating point Python list of the centrality values.
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    percentage: float
        Values between [0, 1]. Default value is 0.10, which means that 10
        percent of the residues will be plotted as VdW spheres.

    Returns
    -------
    Nothing

    """

    numKeyResidues = int(percentage * len(selectedAtoms))
    PML_FILE = open(out_file + '_' + centrality + '.pml', 'w')
    
    PML_FILE.write("load " + out_file + "_" + centrality + ".pdb" + "\n")
    PML_FILE.write("cartoon type = tube\n")
    PML_FILE.write("spectrum b\n")
    PML_FILE.write("set sphere_scale, 0.75\n\n")

    vdw_representation_string = "show spheres, chain {0:s} and resi {1:d} and name ca\n"

    sortedList = np.flip(np.argsort(centralityArray))
    for i in range(0, numKeyResidues):
        # print(centralityArray[sortedList[i]])
        PML_FILE.write(vdw_representation_string.\
            format(selectedAtoms.getChids()[sortedList[i]],
                    selectedAtoms.getResnums()[sortedList[i]]))

    selectedAtoms.setBetas([ i for i in centralityArray])

    writePDB(out_file + '_' + centrality + '.pdb', selectedAtoms)

    PML_FILE.close()

def plotCentralities(centrality, centralityArray, out_file, selectedAtoms, \
                    scalingFactor):
    """
    Plots the centrality values on a 2D graph.

    The centrality values are plotted on a 2D png file.
    If there are at least two chains, the function produces a figure
    for each chain.

    Parameters
    ----------
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness' or 'current_flow_closeness'.
    centralityArray: A numpy data array ?
        It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    ScalingFactor: float
        Sometimes, the values of the centrality arrays are too small.
        The scaling factor multiplies the array to make the values visible in
        the Bfactor colums.

    Return
    ------
    Nothing

    """

    chains = Counter(selectedAtoms.getChids()).keys()

    if len(chains) > 1:
        # Intra-chain
        for chain in chains:
            x = []
            y = []
            dst_file = f"{out_file}_{centrality}_chain{chain}"
            for i in range(0, len(selectedAtoms.getResnums())):
                if selectedAtoms.getChids()[i] == chain:
                    x.append(selectedAtoms[i].getResnum())
                    y.append(centralityArray[i])

            fig, ax = plt.subplots()
            plt.title('Chain ' + chain)
            plt.locator_params(axis='y', nbins=5)

            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            plt.ylabel(centrality.replace('_', ' '), fontsize=20)
            plt.xlabel("Residue Number", fontsize=20)

            # plt.plot(x, centralityArray '.', color='k')
            plt.bar(x, y, color='k')
            plt.tight_layout()
            # plt.show()
            plt.savefig(dst_file + '.png')
            plt.close('all')
    else:
        dst_file = f"{out_file}_{centrality}"
        fig, ax = plt.subplots()
        # plt.locator_params(axis='y', nbins=4)

        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        plt.ylabel(centrality.replace('_', ' '), fontsize=20)
        plt.xlabel("Residue Number", fontsize=20)

        x = selectedAtoms.getResnums()
        # plt.plot(x, centralityArray '.', color='k')
        plt.bar(x, centralityArray, color='k')
        plt.tight_layout()
        # plt.xlim(xmin=0) #xmin will be removed in matplotlib 3.2
        # plt.xlim(left=0)
        # plt.show()
        plt.savefig(dst_file + '.png')
        plt.close('all')


def projectCommunitiesOntoProteinVMD(sortedCommunities, out_file, selectedAtoms):
    """
    Produces VMD output files for visualizing protein communities.

    This function writes a tcl file and a PDB file that can be viewed in
    VMD. Occupancy field of the protein contains the community information.

    The output files can be visualized with VMD (Visual Molecular
    dynamics) program as follows.
    i) Load your pdb file, whether via command line or graphical interface.
    ii) Go to Extensions -> Tk Console and then
    iii) source vmd-output-general.tcl

    Parameters
    ----------
    sortedCommunities: Iterator over tuples of sets of nodes
        It is a tuple of lists. Each list contain a community.
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    Returns
    -------
    Nothing

    """
    VMD_FILE = open(out_file + '_communities.tcl', 'w')

    VMD_FILE.write("mol new " + out_file + '_communities.pdb' + "\n")
    VMD_FILE.write("mol modstyle 0 0 Tube 0.5 25\n")
    VMD_FILE.write("mol modcolor 0 0 Occupancy\n")
    VMD_FILE.write("mol modmaterial 0 0 Glossy\n")

    selectedAtoms.setOccupancies(0)
    i=0
    for item in sortedCommunities:
        #print(item)
        for residx in item:
            selectedAtoms[residx].setOccupancy(i)
        i = i + 1  

    writePDB(out_file + '_communities.pdb', selectedAtoms)

    VMD_FILE.close()

def projectCommunitiesOntoProteinPyMol(sortedCommunities, out_file, selectedAtoms):
    """
    Produces PyMol output files for visualizing protein communities.

    This function writes a pml file and a PDB file that can be viewed in
    PyMol. Occupancy field of the protein contains the community information.

    The output files can be visualized with PyMol program as follows:
    pymol outputfile.pml

    Parameters
    ----------
    sortedCommunities: Iterator over tuples of sets of nodes
        It is a tuple of lists. Each list contain a community.
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    Returns
    -------
    Nothing

    """
    PML_FILE = open(out_file + '_communities.pml', 'w')

    PML_FILE.write("load " + out_file + '_communities.pdb' + "\n")
    PML_FILE.write("cartoon type = tube\n")
    PML_FILE.write("spectrum q\n")
    PML_FILE.write("set sphere_scale, 0.75\n\n")

    selectedAtoms.setOccupancies(0)
    i=0
    for item in sortedCommunities:
        #print(item)
        for residx in item:
            selectedAtoms[residx].setOccupancy(i)
        i = i + 1  

    writePDB(out_file + '_communities.pdb', selectedAtoms)

    PML_FILE.close()

def buildDynamicsNetwork(ccMatrix, distanceMatrix, \
                        valueFilter, distanceFilter,\
                        selectedAtoms):
    """
    This function build a network (graph) from a dynamical correlation
    matrix. 

    The C_ij correlation values are converted to network edges according 
    to -log(abs(C_ij)) (See https://doi.org/10.1073/pnas.0810961106 for 
    details). When only C_ij values are between [0-1], it gives consistent 
    results.

    Parameters
    ----------
    ccMatrix: Numpy matrix
        It is a numpy matrix of typically nDCC, nLMI or Generalized Correlations.
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

def buildSequenceNetwork(ccMatrix, distanceMatrix, \
                        valueFilter, distanceFilter,\
                        selectedAtoms):
    """
    This function build a network (graph) from a dynamical correlation
    matrix. 

    The C_ij correlation values are converted to network edges according 
    to (1.0/abs(C_ij)). It diminishes very fast but it doesn't give negative 
    weigths for values greater than 1.0. 

    Parameters
    ----------
    ccMatrix: Numpy matrix
        It is a numpy matrix of typically DCC, LMI or any other matrix where 
        absolute values correlations are not between zero and one.
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
    seqNetwork = nx.Graph()

    n = selectedAtoms.numAtoms()

    # Add all CA atoms as nodes
    for i in range(n):
        seqNetwork.add_node(i)

    # Add all pairwise interactions greater than the valueFilter as edges.
    # In addition, add only edges which has a distance of lower than the 
    # distance filter
    for i in range(n):
        for j in range(n):
            if fabs(ccMatrix[i][j]) > valueFilter and distanceMatrix[i][j] <= distanceFilter:
                seqNetwork.add_edge(i, j, weight=(1.0/fabs(ccMatrix[i][j])))
 
    return seqNetwork

def calcEigenvector2Centrality(ccMatrix, distanceMatrix, 
                                selectedAtoms, out_file, localityFactor=5.0):
    """
    This function build a network (graph) from a dynamical correlation
    matrix and a distance matrix. 

    The network is build as explained in the 
    https://www.pnas.org/doi/epdf/10.1073/pnas.1810452115

    Parameters
    ----------
    ccMatrix: Numpy matrix
        It is a numpy matrix of typically nDCC, nLMI or Generalized Correlations.
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
        based on contacts and their preservation during the entire simulation.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.

    Returns
    -------
    eigenvectorResult

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
            dynNetwork.add_edge(i, j, weight=(-log(fabs(ccMatrix[i][j])) * np.exp(-distanceMatrix[i][j])/localityFactor) )
            # if fabs(ccMatrix[i][j]) > valueFilter and distanceMatrix[i][j] <= distanceFilter:
            #     dynNetwork.add_edge(i, j, weight=-log(fabs(ccMatrix[i][j])))
            #     # dynNetwork.add_edge(i, j, weight=fabs(correlationArray[i][j]))
    
    eigenvectorResult = nx.eigenvector_centrality_numpy(dynNetwork, weight='weight')

    print("@> Eigenvector2 calculation finished!")

    # open a file for eigenvectors
    eigenvectorFile = open(f"{out_file}_eigenvector2_locality_factor{localityFactor:.2f}.dat", "w")

    for i in range(n):
        #    print(str(i)+" "+(str(graph.closeness(i, weight='weight'))))
        eigenvectorFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                eigenvectorResult[i],
                                                                selectedAtoms[i].getChid()))
    eigenvectorFile.close()
    centralityArrayScaled = autoScaleCentralities(list(eigenvectorResult.values()),
                                    selectedAtoms, 'local')
    projectCentralitiesOntoProteinVMD('eigenvector2',
                                        centralityArrayScaled,
                                        out_file,
                                        selectedAtoms, scalingFactor=1)
    projectCentralitiesOntoProteinPyMol('eigenvector2',
                                        centralityArrayScaled,
                                        out_file,
                                        selectedAtoms, scalingFactor=1)
    plotCentralities('eigenvector2', centralityArrayScaled, out_file,
                                        selectedAtoms, scalingFactor=1)

    return eigenvectorResult

def centralityAnalysis(graph, valueFilter, distanceFilter, \
                        out_file, centrality, selectedAtoms):
    """
    This function calculates various network (graph) centralities of a protein.

    This function calculates some network centrality measures such as
    degree, betweenness, closeness, current flow betweenness and eigenvector.
    This function needs Python 3.6 or later to maintain dictionary order.!!!

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
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness', 'current_flow_closeness', 'eigenvector'
        or 'community'.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.

    Returns
    -------
    Nothing

    """
    # dynNetwork = buildDynamicsNetwork(ccMatrix, distanceMatrix, \
    #                                 valueFilter, distanceFilter,\
    #                                 selectedAtoms)

    n = selectedAtoms.numAtoms()
    ########################## Calculate degrees of all nodes
    if centrality == 'degree':
        degreeResult = graph.degree(weight='weight')
        degreeResultList = []
        for i in range(0, len(degreeResult)):
            degreeResultList.append(degreeResult[i])
        # print(degreeResultList)
        print("@> Degree calculation finished!")

        # open a file for degree
        degreeFile = open(f"{out_file}_degree_value_filter{valueFilter:.2f}.dat", "w")
        for i in range(n):
            #    print(str(i)+" "+(str(graph.degree(i, weight='weight'))))
            degreeFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                              degreeResult[i],
                                                              selectedAtoms[i].getChid()))
        degreeFile.close()
        projectCentralitiesOntoProteinVMD(centrality,
                                          degreeResultList,
                                          out_file,
                                          selectedAtoms,
                                          scalingFactor=1)
        projectCentralitiesOntoProteinPyMol(centrality,
                                          degreeResultList,
                                          out_file,
                                          selectedAtoms,
                                          scalingFactor=1)
        plotCentralities(centrality,
                         degreeResultList,
                         out_file,
                         selectedAtoms,
                         scalingFactor=1)


    ########################## Calculate betweenness
    elif centrality == 'betweenness':
        betweennessResult = nx.betweenness_centrality(graph, k=None,
                                                      normalized=True, weight='weight',
                                                      endpoints=False, seed=None)
        print("@> Betweenness calculation finished!")
        # print(betweennessResult)

        # open a file for betweenness
        betweennessFile = open(f"{out_file}_betweenness_value_filter{valueFilter:.2f}.dat", "w")

        for i in range(n):
            betweennessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".\
                format(selectedAtoms[i].getResnum(),
                        betweennessResult[i],
                        selectedAtoms[i].getChid()))
        betweennessFile.close()
        # print(list(betweennessResult.values()))
        # print(len(list(betweennessResult.values())))
        projectCentralitiesOntoProteinVMD(centrality,
                                          list(betweennessResult.values()),
                                          out_file, selectedAtoms, scalingFactor=1000)
        projectCentralitiesOntoProteinPyMol(centrality,
                                          list(betweennessResult.values()),
                                          out_file, selectedAtoms, scalingFactor=1000)
        plotCentralities(centrality,
                         list(betweennessResult.values()),
                         out_file, selectedAtoms, scalingFactor=1)

    ########################## Calculate closeness
    elif centrality == 'closeness':
        closenessResult = nx.closeness_centrality(graph, u=None, distance='weight')
        print("@> Closeness calculation finished!")

        # open a file for closeness
        closenessFile = open(out_file + "_closeness_value_filter" + "{:.2f}".format(valueFilter) + '.dat', "w")

        for i in range(n):
            #    print(str(i)+" "+(str(graph.closeness(i, weight='weight'))))
            closenessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                 closenessResult[i],
                                                                 selectedAtoms[i].getChid()))
        closenessFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(closenessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1)
        projectCentralitiesOntoProteinPyMol(centrality,
                                          list(closenessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1)
        plotCentralities(centrality,
                         list(closenessResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)


    ########################## Calculate current_flow_betweenness
    elif centrality == 'current_flow_betweenness':
        current_flow_betweennessResult = nx.current_flow_betweenness_centrality(graph, normalized=True,
                                                                                weight='weight')

        print("@> Current flow betweenness calculation finished!")

        # open a file for current_flow betweenness
        current_flow_betweennessFile = open(f"{out_file}_current_flow_betweenness_value_filter{valueFilter:.2f}.dat",
                                            "w")

        for i in range(n):
            #    print(str(i)+" "+(str(graph.betweenness(i, weight='weight'))))
            current_flow_betweennessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".\
                format(selectedAtoms[i].getResnum(),
                        current_flow_betweennessResult[i],
                        selectedAtoms[i].getChid()))
        current_flow_betweennessFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(current_flow_betweennessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1000)
        projectCentralitiesOntoProteinPyMol(centrality,
                                          list(current_flow_betweennessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1000)
        plotCentralities(centrality,
                         list(current_flow_betweennessResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)
    ########################## Calculate closeness
    elif centrality == 'current_flow_closeness':
        current_flow_closenessResult = nx.current_flow_closeness_centrality(graph, weight='weight')

        print("@> Current flow closeness calculation finished!")

        # open a file for current_flow closeness
        current_flow_closenessFile = open(f"{out_file}_current_flow_closeness_value_filter{valueFilter:.2f}.dat",
                                          "w")

        for i in range(n):
            #    print(str(i)+" "+(str(graph.closeness(i, weight='weight'))))
            current_flow_closenessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".\
                format(selectedAtoms[i].getResnum(),
                        current_flow_closenessResult[i],
                        selectedAtoms[i].getChid()))
        current_flow_closenessFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(current_flow_closenessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1000)
        projectCentralitiesOntoProteinPyMol(centrality,
                                          list(current_flow_closenessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1000)
        plotCentralities(centrality,
                         list(current_flow_closenessResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)
    ########################## Calculate eigenvector centrality
    elif centrality == 'eigenvector':
        eigenvectorResult = nx.eigenvector_centrality_numpy(graph, weight='weight')
        print("@> Eigenvector calculation finished!")

        # open a file for eigenvectors
        eigenvectorFile = open(f"{out_file}_eigenvector_value_filter{valueFilter:.2f}.dat", "w")

        for i in range(n):
            #    print(str(i)+" "+(str(graph.closeness(i, weight='weight'))))
            eigenvectorFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                   eigenvectorResult[i],
                                                                   selectedAtoms[i].getChid()))
        eigenvectorFile.close()

        # centralityArrayScaled = autoScaleCentralities(centrality, \
        #                                     list(eigenvectorResult.values()),
        #                                     selectedAtoms, 'global')
        projectCentralitiesOntoProteinVMD(centrality,
                                            list(eigenvectorResult.values()),
                                            out_file,
                                            selectedAtoms, scalingFactor=1)
        projectCentralitiesOntoProteinPyMol(centrality,
                                            list(eigenvectorResult.values()),
                                            out_file,
                                            selectedAtoms, scalingFactor=1)
        plotCentralities(centrality, list(eigenvectorResult.values()), out_file,
                                            selectedAtoms, scalingFactor=1)
    ########################## Calculate communities with Girvan-Newman
    elif centrality == 'community':
        from networkx.algorithms import community
        communities = community.girvan_newman(graph)
        sortedCommunities = tuple(sorted(c) for c in next(communities))
        print("@> There are " + str(len(sortedCommunities)) + \
                " communities in your structure.")
        projectCommunitiesOntoProteinVMD(sortedCommunities, out_file, selectedAtoms)
        projectCommunitiesOntoProteinPyMol(sortedCommunities, out_file, selectedAtoms)
        print("@> Community calculation finished!")

    else:
        print("ERROR: Unknown centrality selected! It can only be")
        print("       'degree', 'betweenness', 'closeness',")
        print("       'current_flow_betweenness', 'current_flow_closeness'")
        print("       or 'eigenvector!'")
        sys.exit(-1)

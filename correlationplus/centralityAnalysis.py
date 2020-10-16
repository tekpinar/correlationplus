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

import sys
from math import fabs, log
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from prody import writePDB
import networkx as nx


def projectCentralitiesOntoProteinVMD(centrality, centralityArray, out_file, selectedAtoms, scalingFactor):
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
    ii) Go to Extemsions -> Tk Console and then
    iii) source vmd-output-general.tcl
    It can take some to load the general script.

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


def plotCentralities(centrality, centralityArray, out_file, selectedAtoms, scalingFactor):
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
            dst_file = out_file + '_' + centrality + '_chain' + chain
            for i in range(0, len(selectedAtoms.getResnums())):
                if selectedAtoms.getChids()[i] == chain:
                    x.append(selectedAtoms[i].getResnum())
                    y.append(centralityArray[i])

            plt.subplots()
            plt.locator_params(axis='y', nbins=5)

            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.ylabel(centrality.replace('_', ' '), fontsize=20)
            plt.xlabel("Residue Number", fontsize=20)

            # plt.plot(x, centralityArray '.', color='k')
            plt.bar(x, y, color='k')
            plt.tight_layout()
            # plt.show()
            plt.savefig(dst_file + '.png')
            plt.close('all')
    else:
        dst_file = out_file + '_' + centrality
        plt.subplots()
        # plt.locator_params(axis='y', nbins=4)

        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
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


def centralityAnalysis(ccMatrix, distanceMatrix, valueFilter, distanceFilter, out_file, centrality, selectedAtoms):
    """
    This function calculates various network (graph) centralities of a protein.

    This function calculates some network centrality measures such as
        -degree
        -betweenness
        -closeness
        -current flow betweenness
        -eigenvector.
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
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be
        extended.
    centrality: string
        It can have 'degree', 'betweenness', 'closeness',
        'current_flow_betweenness' or 'current_flow_closeness'.
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.

    Returns
    -------
    Nothing
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
            if((fabs(ccMatrix[i][j])>valueFilter) and (distanceMatrix[i][j]<=distanceFilter)):
                dynNetwork.add_edge(i, j, weight=-log(fabs(ccMatrix[i][j])))
                # dynNetwork.add_edge(i, j, weight=fabs(correlationArray[i][j]))

    ##########################Calculate degrees of all nodes
    if centrality == 'degree':
        degreeResult = dynNetwork.degree(weight='weight')
        degreeResultList = []
        for i in range(0, len(degreeResult)):
            degreeResultList.append(degreeResult[i])
        # print(degreeResultList)
        print("@> Degree calculation finished!")

        # open a file for degree
        degreeFile = open(out_file + "_degree_value_filter" + "{:.2f}".format(valueFilter) + '.dat', "w")
        for i in range(n):
            #    print(str(i)+" "+(str(dynNetwork.degree(i, weight='weight'))))
            degreeFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                              degreeResult[i],
                                                              selectedAtoms[i].getChid()))
        degreeFile.close()
        projectCentralitiesOntoProteinVMD(centrality,
                                          degreeResultList,
                                          out_file,
                                          selectedAtoms,
                                          scalingFactor=1)
        plotCentralities(centrality,
                         degreeResultList,
                         out_file,
                         selectedAtoms,
                         scalingFactor=1)


    ##########################Calculate betweenness
    elif centrality == 'betweenness':
        betweennessResult = nx.betweenness_centrality(dynNetwork, k=None,
                                                      normalized=True, weight='weight',
                                                      endpoints=False, seed=None)
        print("@> Betweenness calculation finished!")
        # print(betweennessResult)

        # open a file for betweenness
        betweennessFile = open(out_file + "_betweenness_value_filter" + "{:.2f}".format(valueFilter) + '.dat', "w")

        for i in range(n):
            #    print(str(i)+" "+(str(dynNetwork.betweenness(i, weight='weight'))))
            betweennessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                   betweennessResult[i],
                                                                   selectedAtoms[i].getChid()))
        betweennessFile.close()
        # print(list(betweennessResult.values()))
        # print(len(list(betweennessResult.values())))
        projectCentralitiesOntoProteinVMD(centrality,
                                          list(betweennessResult.values()),
                                          out_file, selectedAtoms, scalingFactor=1000)
        plotCentralities(centrality,
                         list(betweennessResult.values()),
                         out_file, selectedAtoms, scalingFactor=1)

    ##########################Calculate closeness
    elif centrality == 'closeness':
        closenessResult = nx.closeness_centrality(dynNetwork, u=None, distance='weight')
        print("@> Closeness calculation finished!")

        # open a file for closeness
        closenessFile = open(out_file + "_closeness_value_filter" + "{:.2f}".format(valueFilter) + '.dat', "w")

        for i in range(n):
            #    print(str(i)+" "+(str(dynNetwork.closeness(i, weight='weight'))))
            closenessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                 closenessResult[i],
                                                                 selectedAtoms[i].getChid()))
        closenessFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(closenessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1)
        plotCentralities(centrality,
                         list(closenessResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)


    ##########################Calculate current_flow_betweenness
    elif centrality == 'current_flow_betweenness':
        current_flow_betweennessResult = nx.current_flow_betweenness_centrality(dynNetwork, normalized=True,
                                                                                weight='weight')

        print("@> Current flow betweenness calculation finished!")

        # open a file for current_flow betweenness
        current_flow_betweennessFile = open(f"{out_file}_current_flow_betweenness_value_filter{valueFilter:.2f}.dat",
                                            "w")

        for i in range(n):
            #    print(str(i)+" "+(str(dynNetwork.betweenness(i, weight='weight'))))
            current_flow_betweennessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                                current_flow_betweennessResult[i],
                                                                                selectedAtoms[i].getChid()))
        current_flow_betweennessFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(current_flow_betweennessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1000)
        plotCentralities(centrality,
                         list(current_flow_betweennessResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)
    ##########################Calculate closeness
    elif centrality == 'current_flow_closeness':
        current_flow_closenessResult = nx.current_flow_closeness_centrality(dynNetwork, weight='weight')

        print("@> Current flow closeness calculation finished!")

        # open a file for current_flow closeness
        current_flow_closenessFile = open(f"{out_file}_current_flow_closeness_value_filter{valueFilter:.2f}.dat",
                                          "w")

        for i in range(n):
            #    print(str(i)+" "+(str(dynNetwork.closeness(i, weight='weight'))))
            current_flow_closenessFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                              current_flow_closenessResult[i],
                                                                              selectedAtoms[i].getChid()))
        current_flow_closenessFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(current_flow_closenessResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1000)
        plotCentralities(centrality,
                         list(current_flow_closenessResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)
    ##########################Calculate eigenvector centrality
    elif centrality == 'eigenvector':
        eigenvectorResult = nx.eigenvector_centrality_numpy(dynNetwork, weight='weight')
        print("@> Eigenvector calculation finished!")

        # open a file for closeness
        eigenvectorFile = open(out_file + "_eigenvector_value_filter" + "{:.2f}".format(valueFilter) + '.dat', "w")

        for i in range(n):
            #    print(str(i)+" "+(str(dynNetwork.closeness(i, weight='weight'))))
            eigenvectorFile.write("{0:d}\t{1:.6f}\t{2:s}\n".format(selectedAtoms[i].getResnum(),
                                                                   eigenvectorResult[i],
                                                                   selectedAtoms[i].getChid()))
        eigenvectorFile.close()

        projectCentralitiesOntoProteinVMD(centrality,
                                          list(eigenvectorResult.values()),
                                          out_file,
                                          selectedAtoms, scalingFactor=1)
        plotCentralities(centrality,
                         list(eigenvectorResult.values()),
                         out_file,
                         selectedAtoms, scalingFactor=1)
    else:
        print("ERROR: Unknown centrality selected! It can only be")
        print("       'degree', 'betweenness', 'closeness',")
        print("       'current_flow_betweenness', 'current_flow_closeness'")
        print("       or 'eigenvector!'")
        sys.exit(-1)

##############################################################################
# correlationPlus - Python module to plot dynamical correlations maps         #
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

from collections import Counter, OrderedDict

import matplotlib.pyplot as plt
import numpy as np

from .mapAnalysis import cmap_discretize

def overallUniformDifferenceMap(ccMatrix1, ccMatrix2,
                                minColorBarLimit, maxColorBarLimit, 
                                out_file, title, selectedAtoms):
    """
        Plots the difference map between correlation maps for the entire structure.
        Sizes of ccMatrix1 and ccMatrix2 are identical. Only one atom set is
        sufficient to plot the difference map.
        
    Parameters
    ----------
    ccMatrix1: A numpy square matrix of floats
        The first correlation matrix.
    ccMatrix2: A numpy square matrix of floats
        The second correlation matrix.
    minColorBarLimit: signed int
        If nDCC maps -2. If absndcc or lmi, -1.  
    maxColorBarLimit: unsigned int
        If nDCC maps 2. If absndcc or lmi, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.

    Returns
    -------
    Nothing
    """

    # selectedAtoms = parsePDB(pdb_file, subset='ca')
    # chainLengths = Counter(selectedAtoms.getChids()).values()
    diffMap = np.subtract(ccMatrix1, ccMatrix2)

    n = len(ccMatrix1)
    ##########################################################################
    # Set residue interface definitions
    fig = plt.figure()
    fig.set_size_inches(8.0, 5.5, forward=True)
    plt.rcParams['font.size'] = 16
    ax = fig.add_subplot(1, 1, 1)

    plt.xlabel('Residue indices')
    plt.ylabel('Residue indices')
    plt.title(title, y=1.08)
    plt.grid(color='w', linestyle='--', linewidth=1)

    # print(selectedAtoms.getChids())

    myList = list(Counter(selectedAtoms.getChids()).keys())

    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtoms.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    for i in Counter(selectedAtoms.getChids()).values():
        tempVal = tempVal + i
        selection_reorder.append(tempVal)
        selection_tick_labels.append(str(selectedAtoms.getResnums()[tempVal - 1]))

    # print(selection_reorder)

    major_nums = []
    major_labels = []
    major_nums.extend(selection_reorder)
    major_labels.extend(selection_tick_labels)

    ##########################################################################
    # Set plotting parameters

    # plt.rcParams['axes.titlepad'] = 20
    ax.autoscale(False)
    ax.set_aspect('equal')

    # ax.set_xticks(major_nums, major_labels, rotation=45, minor=False)
    plt.xticks(major_nums, major_labels, size=12, rotation=45)
    plt.yticks(major_nums, major_labels, size=12)

    # ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
    # ax.yaxis.set_tick_params(width=2, length=5)

    plt.axis([0, n, 0, n])
    ax.tick_params(which='major', width=2, length=5)
    ax.tick_params(which='minor', width=1, length=3)

    # print("Min value")
    # print("Row Index of min value")
    # print(np.argmin(ccMatrix1_sub, 0))

    # print("Column Index of min value")
    # print(np.argmin(ccMatrix1_sub, 1))

    # Set colorbar features here!
    jet = plt.get_cmap('jet')
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(diffMap), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position = fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar = plt.colorbar(cax=position)

    cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')
        t.set_x(4.0)

    # Add chain borders to the plot
    for i in range(len(selection_reorder) - 1):
        beginningPoint = selection_reorder[i] / selection_reorder[-1]
        endingPoint = selection_reorder[i + 1] / selection_reorder[-1]
        middlePoint = (float(beginningPoint) + float(endingPoint)) / 2.0
        if i % 2 == 0:
            # x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction',
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth=2., arrowstyle="-", color='black'))

            # y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction',
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth=2., arrowstyle="-", color='black'))

        elif i % 2 == 1:
            # x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction',
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth=2., arrowstyle="-", color='gray'))

            # y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction',
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth=2., arrowstyle="-", color='gray'))

        ax.annotate(myList[i], xy=(0, 1.04), xycoords='axes fraction', xytext=(middlePoint - 0.015, 1.04),
                    size=14, color='black')
        ax.annotate(myList[i], xy=(1.05, 0), xycoords='axes fraction', xytext=(1.05, middlePoint - 0.015),
                    rotation=90, size=14, color='black')
        # print(middlePoint)

    # plt.tight_layout()
    plt.savefig(out_file + '-overall-difference.png', bbox_inches='tight', dpi=200)
    # plt.show()


def overallNonUniformDifferenceMap(ccMatrix1, ccMatrix2, minColorBarLimit,
                                   maxColorBarLimit, out_file, title,
                                   selectedAtomSet1, selectedAtomSet2):

    """
        Plots the difference map between correlation maps for the entire structure.
        Sizes of ccMatrix1 and ccMatrix2 are not identical. A mapping for matching
        residues is performed before difference map plotting.

    Parameters
    ----------
    ccMatrix1: A numpy square matrix of floats
        The first correlation matrix.
    ccMatrix2: A numpy square matrix of floats
        The second correlation matrix.
    minColorBarLimit: signed int
        If nDCC maps -2. If absndcc or lmi, -1.  
    maxColorBarLimit: unsigned int
        If nDCC maps 2. If absndcc or lmi, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtomSet1: prody object
        A list of -typically CA- atoms selected from the first PDB file.
    selectedAtomSet2: prody object
        A list of -typically CA- atoms selected from the second PDB file.

    Returns
    -------
    Nothing
    """

    # selectedAtoms = parsePDB(pdb_file, subset='ca')
    # chainLengths = Counter(selectedAtoms.getChids()).values()

    commonCoreDictionary = findCommonCorePDB(selectedAtomSet1, selectedAtomSet2)
    diffMap = np.zeros((len(commonCoreDictionary), len(commonCoreDictionary)), dtype=float)

    items = list(commonCoreDictionary.items())

    for i in range(0, len(commonCoreDictionary)):
        key1, value1 = items[i]
        for j in range(i, len(commonCoreDictionary)):
            key2, value2 = items[j]
            diffMap[i][j] = ccMatrix1[key1][key2] - ccMatrix2[value1][value2]
            diffMap[j][i] = diffMap[i][j]
            # diffMap = np.subtract(ccMatrix1, ccMatrix2)

    n = len(commonCoreDictionary)
    ##########################################################################
    # Set residue interface definitions
    fig = plt.figure()
    fig.set_size_inches(8.0, 5.5, forward=True)
    plt.rcParams['font.size'] = 16
    ax = fig.add_subplot(1, 1, 1)

    plt.xlabel('Residue indices')
    plt.ylabel('Residue indices')
    plt.title(title, y=1.08)
    plt.grid(color='w', linestyle='--', linewidth=1)

    ##########################################################################
    # Set xtics, ytic, xlabels, y labels etc.
    # This part may have to be rewritten!
    myList = list(Counter(selectedAtomSet1.getChids()).keys())

    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtomSet1.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    for i in Counter(selectedAtomSet1.getChids()).values():
        tempVal = tempVal + i
        selection_reorder.append(tempVal)
        selection_tick_labels.append(str(selectedAtomSet1.getResnums()[tempVal - 1]))

    # print(selection_reorder)

    major_nums = []
    major_labels = []
    major_nums.extend(selection_reorder)
    major_labels.extend(selection_tick_labels)

    ##########################################################################
    # Set plotting parameters

    # plt.rcParams['axes.titlepad'] = 20
    ax.autoscale(False)
    ax.set_aspect('equal')

    # ax.set_xticks(major_nums, major_labels, rotation=45, minor=False)
    plt.xticks(major_nums, major_labels, size=12, rotation=45)
    plt.yticks(major_nums, major_labels, size=12)

    # ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
    # ax.yaxis.set_tick_params(width=2, length=5)

    plt.axis([0, n, 0, n])
    ax.tick_params(which='major', width=2, length=5)
    ax.tick_params(which='minor', width=1, length=3)

    # print("Min value")
    # print("Row Index of min value")
    # print(np.argmin(ccMatrix1_sub, 0))

    # print("Column Index of min value")
    # print(np.argmin(ccMatrix1_sub, 1))

    # Set colorbar features here!
    jet = plt.get_cmap('jet')
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(diffMap), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position = fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar = plt.colorbar(cax=position)

    cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')
        t.set_x(4.0)
    ###########################################################################
    # Add chain borders to the plot
    for i in range(len(selection_reorder) - 1):
        beginningPoint = selection_reorder[i] / selection_reorder[-1]
        endingPoint = selection_reorder[i + 1] / selection_reorder[-1]
        middlePoint = (float(beginningPoint) + float(endingPoint)) / 2.0
        if (i % 2 == 0):
            # x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction',
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth=2., arrowstyle="-", color='black'))

            # y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction',
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth=2., arrowstyle="-", color='black'))

        elif i % 2 == 1:
            # x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction',
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth=2., arrowstyle="-", color='gray'))

            # y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction',
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth=2., arrowstyle="-", color='gray'))

        ax.annotate(myList[i], xy=(0, 1.04),
                    xycoords='axes fraction',
                    xytext=(middlePoint - 0.015, 1.04),
                    size=14, color='black')
        ax.annotate(myList[i], xy=(1.05, 0),
                    xycoords='axes fraction',
                    xytext=(1.05, middlePoint - 0.015),
                    rotation=90,
                    size=14, color='black')
        # print(middlePoint)
    ###########################################################################
    # plt.tight_layout()
    plt.savefig(out_file + '-overall-difference.png', bbox_inches='tight', dpi=200)
    # plt.show()


def findCommonCorePDB(selectedAtomSet1, selectedAtomSet2):
    """
        Finds a common set of residues between different conformations of a 
        protein.
        
        This function assumes that two structures are obtained
        exactly from the same species and there is not any mutation in any one
        of them. This can happen when two conformations of a protein are obtained
        with different number of missing atoms.
        Under these conditions, we are trying to match two structures and find
        the common core of two structures that have different number of CA atoms.
        We will use this information to subtract two correlation maps!
        
    Parameters
    ----------
    ccMatrix1: A numpy square matrix of floats
        The first correlation matrix.
    ccMatrix2: A numpy square matrix of floats
        The second correlation matrix.
    minColorBarLimit: signed int
        If nDCC maps -2. If absndcc or lmi, -1.  
    maxColorBarLimit: unsigned int
        If nDCC maps 2. If absndcc or lmi, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtomSet1: prody object
        A list of -typically CA- atoms selected from the first PDB file.
    selectedAtomSet2: prody object
        A list of -typically CA- atoms selected from the second PDB file.

    Returns
    -------
    commonCoreDictionary: dictionary
        A dictionary of residues in conformation A and corresponding residue in 
        conformation B. 
    
    
    """
    print("@> Calculating a common core for two structures.")
    lengthSet1 = len(selectedAtomSet1)
    lengthSet2 = len(selectedAtomSet2)

    # print(lengthSet1)
    # print(lengthSet2)

    # I made it an ordered dictionary so that the it will not shuffle
    # the items in the dictionary!
    commonCoreDictionary = OrderedDict()
    # if (lengthSet1 > lengthSet2):
    for i in range(0, lengthSet1):
        for j in range(0, lengthSet2):
            if (selectedAtomSet1.getResnums()[i] == selectedAtomSet2.getResnums()[j]) and \
                    (selectedAtomSet1.getResnames()[i] == selectedAtomSet2.getResnames()[j]) and \
                    (selectedAtomSet1.getChids()[i] == selectedAtomSet2.getChids()[j]):
                commonCoreDictionary[i] = j
    # print(commonCoreDictionary)
    return commonCoreDictionary

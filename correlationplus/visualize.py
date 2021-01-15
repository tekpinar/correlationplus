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

from collections import Counter, OrderedDict
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
from prody import buildDistMatrix



def cmap_discretize(cmap, N):
    """
        Creates a discrete colormap from the continuous colormap cmap.

    Parameters
    ----------

    cmap: colormap instance, eg. cm.jet.
    N: number of colors.

    Returns
    -------
    cmap: A discrete color map.     
    
    Example
    -------
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1, ki], colors_rgba[i, ki]) for i in range(N + 1)]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)


def convertLMIdata2Matrix(inp_file, writeAllOutput: bool):
    """
        This function parses LMI matrix and returns a numpy array. If the 
        It can handle both full matrix format or g_correlation format.
    Parameters
    ----------

    inp_file: string
        LMI file to read.
 
    writeAllOutput: bool
        If True, an output file for the LMI values will be written in matrix 
        format. The matrix does not contain residue names etc.  

    Returns
    -------
    cc: A numpy array of float value arrays.
        LMI values in matrix format.  

    """
    data_file = open(inp_file, 'r')
    allLines = data_file.readlines()
    data_file.close()
    if(("x" in allLines[0]) and ("[" in allLines[0])):
        # This is a g_correlation file.
        data_list = []
        for line in allLines:
            words = line.split()
            for i in words:
                    data_list.append(i)

        data_list = data_list[4:-1]
        n = int(np.sqrt(len(data_list)))
        data_array = np.array(data_list, dtype=float)

        cc = np.reshape(data_array, (n, n))

        if writeAllOutput:
            np.savetxt(inp_file[:-4] + "_modif.dat", cc, fmt='%.6f')
    else:
        # The data is assumed to be in full matrix format.
        # Please note that correlationplus can not handle upper or lower
        # triangle matrix format for now.
        cc = np.loadtxt(inp_file, dtype=float)

    return cc


def filterCorrelationMapByDistance(ccMatrix, out_file, title,
                                   selectedAtoms, distanceValue,
                                   absoluteValues: bool, writeAllOutput: bool):

    
    """
        If residues are closer to each other than a certain distance
    (distanceValue), make these correlations zero. This filter can be useful
    to get rid of high short distance correlations for visualization purposes.
    This function returns a filtered ccMatrix.

    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    distanceValue: float
        A distance value in Angstrom unit. For example, it is good to remove 
        high correlations for residues within less 5.0 Angstrom distance 
        to have a clear visualization.
    absoluteValues: bool
        If True, an absolute values file will be written. 
    writeAllOutput: bool
        If True, an output file for distances and correlation values can
        be written. This can be useful to see their distribution as well as 
        individual values. 

    Returns
    -------
    Nothing
    """
    print("@> Filtering correlations lower than " + str(distanceValue) + " Angstrom.")
    # Calculate distance matrix
    dist_matrix = buildDistMatrix(selectedAtoms)

    # print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    # print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    for i in range(0, len(ccMatrix)):
        for j in range(i + 1, len(ccMatrix)):
            if dist_matrix[i][j] < distanceValue:
                ccMatrix[i][j] = 0.0
                ccMatrix[j][i] = 0.0

    if absoluteValues:
        dst_file = out_file + '-absolute-correlation-filtered'

    else:
        dst_file = out_file + '-correlation-filtered'

    # Write output
    # Writing the output is very important for further analyses such as
    # inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    if (writeAllOutput):
        DATA_FILE = open(dst_file + 'filtered.dat', 'w')
        for i in range(0, len(ccMatrix)):
            for j in range(i + 1, len(ccMatrix)):
                DATA_FILE.write("{0:d}\t{1:s}\t{2:d}\t{3:s}\t{4:.3f}\t{5:.3f}\n".format(selectedAtoms.getResnums()[i],
                                                                                        selectedAtoms.getChids()[i],
                                                                                        selectedAtoms.getResnums()[j],
                                                                                        selectedAtoms.getChids()[j],
                                                                                        dist_matrix[i][j],
                                                                                        ccMatrix[i][j]))
        DATA_FILE.close()
    return ccMatrix


def overallCorrelationMap(ccMatrix,
                          minColorBarLimit, maxColorBarLimit,
                          out_file, title, selectedAtoms):
    """
        Plots nDCC maps for the whole structure.

    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
    minColorBarLimit: signed int
        Mostly, -1 or 0.
    maxColorBarLimit: unsigned int
        Mostly, 1.
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

    n = len(ccMatrix)
    ##########################################################################
    # Set residue interface definitions
    fig = plt.figure()
    fig.set_size_inches(8.0, 5.5, forward=True)
    plt.rcParams['font.size'] = 16
    ax = fig.add_subplot(1, 1, 1)

    plt.xlabel('Residue indices')
    plt.ylabel('Residue indices')
    plt.title(title, y=1.08)

    # print(selectedAtoms.getChids())

    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtoms.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    myList = list(Counter(selectedAtoms.getChids()).keys())

    major_nums = []
    major_labels = []

    if len(myList) == 1:
        realTicsList = np.linspace(0, len(selectedAtoms.getResnums()) - 1, 6, dtype=int)
        major_nums = realTicsList

        for item in (realTicsList):
            # print(selectedAtoms.getResnums()[item])
            major_labels.append(str(selectedAtoms.getResnums()[item]))
    elif (len(myList) > 1):
        for i in Counter(selectedAtoms.getChids()).values():
            tempVal = tempVal + i
            selection_reorder.append(tempVal)
            selection_tick_labels.append(str(selectedAtoms.getResnums()[tempVal - 1]))
        major_nums.extend(selection_reorder)
        major_labels.extend(selection_tick_labels)
    else:
        print("Warning: Unknown chain ID!")
    # print(selection_reorder)

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
    # print(np.argmin(ccMatrix_sub, 0))

    # print("Column Index of min value")
    # print(np.argmin(ccMatrix_sub, 1))

    # Set colorbar features here!
    jet = plt.get_cmap('jet')
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(ccMatrix), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position = fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar = plt.colorbar(cax=position)

    cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')
        t.set_x(4.0)

    if len(myList) > 1:
        # Add chain borders to the plot
        for i in range(len(selection_reorder) - 1):
            beginningPoint = selection_reorder[i] / selection_reorder[-1]
            endingPoint = selection_reorder[i + 1] / selection_reorder[-1]
            middlePoint = (float(beginningPoint) + float(endingPoint)) / 2.0
            if (i % 2 == 0):
                # x axis
                ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction',
                            xytext=(endingPoint, 1.03),
                            arrowprops=dict(linewidth=2., arrowstyle="-", color='black'))

                # y axis
                ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction',
                            xytext=(1.04, endingPoint),
                            arrowprops=dict(linewidth=2., arrowstyle="-", color='black'))

            elif i % 2 == 1:
                # x axis
                ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction',
                            xytext=(endingPoint, 1.03), arrowprops=dict(linewidth=2., arrowstyle="-", color='gray'))

                # y axis
                ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction',
                            xytext=(1.04, endingPoint), arrowprops=dict(linewidth=2., arrowstyle="-", color='gray'))

            ax.annotate(myList[i], xy=(0, 1.04), xycoords='axes fraction',
                        xytext=(middlePoint - 0.015, 1.04), size=14, color='black')
            ax.annotate(myList[i], xy=(1.05, 0), xycoords='axes fraction',
                        xytext=(1.05, middlePoint - 0.015), rotation=90, size=14, color='black')
            # print(middlePoint)

    # plt.tight_layout()
    plt.savefig(out_file + '-overall.png', bbox_inches='tight', dpi=200)
    # plt.show()


def intraChainCorrelationMaps(ccMatrix,
                              minColorBarLimit, maxColorBarLimit,
                              out_file, title, selectedAtoms, saveMatrix):
    """
        Plot intra-chain correlations to different png files,
        if there are at least two chains!

    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    minColorBarLimit: signed int
        Mostly, -1 or 0.
    maxColorBarLimit: unsigned int
        Mostly, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    saveMatrix: bool
        If True, an output file for the correlations will be written
        be written. 

    Returns
    -------
    Nothing
    """

    myList = list(Counter(selectedAtoms.getChids()).keys())
    # print(myList)

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

    for j in range(len(myList)):
        # Set labels
        major_nums = []
        major_labels = []
        realTicsList = np.linspace(selection_reorder[j], selection_reorder[j + 1], 6, dtype=int)
        ticsList = np.linspace(selection_reorder[j] - selection_reorder[j],
                               selection_reorder[j + 1] - selection_reorder[j], 6, dtype=int)
        # selection_reorder[j+1]-
        major_nums = ticsList
        # major_nums = major_nums[0:-1]
        # print(realTicsList)
        realTicsList[-1] = realTicsList[-1] - 1
        for item in realTicsList:
            # print(selectedAtoms.getResnums()[item])
            major_labels.append(str(selectedAtoms.getResnums()[item]))

        # print(major_nums)
        # print(major_labels)
        ##########################################################################
        # Set plotting parameters
        ##########################################################################
        # Set residue interface definitions
        fig = plt.figure()
        fig.set_size_inches(8.0, 5.5, forward=True)
        plt.rcParams['font.size'] = 16
        ax = fig.add_subplot(1, 1, 1)

        plt.xlabel('Residue indices')
        plt.ylabel('Residue indices')
        plt.title('Chain ' + myList[j], y=1.08)

        # plt.rcParams['axes.titlepad'] = 20
        ax.autoscale(False)
        ax.set_aspect('equal')
        # for item in ax.get_xticks():
        #    print(item)

        # ax.set_xticks(major_nums, major_labels)
        plt.xticks(major_nums, major_labels, size=12)
        plt.yticks(major_nums, major_labels, size=12)

        # ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
        # ax.yaxis.set_tick_params(width=2, length=5)

        ax.tick_params(which='major', width=2, length=5)
        ax.tick_params(which='minor', width=1, length=3)

        # print("Min value")
        # print("Row Index of min value")
        # print(np.argmin(ccMatrix_sub, 0))

        # print("Column Index of min value")
        # print(np.argmin(ccMatrix_sub, 1))

        # Set colorbar features here!
        jet = plt.get_cmap('jet')
        djet = cmap_discretize(jet, 8)
        plt.axis([0, (selection_reorder[j + 1] - selection_reorder[j]), 0,
                  (selection_reorder[j + 1] - selection_reorder[j])])
        sub_nDCC_matrix = ccMatrix[selection_reorder[j]: selection_reorder[j + 1],
                          selection_reorder[j]: selection_reorder[j + 1]]
        plt.imshow(np.matrix(sub_nDCC_matrix), cmap=djet)
        plt.clim(minColorBarLimit, maxColorBarLimit)

        position = fig.add_axes([0.85, 0.15, 0.03, 0.70])

        cbar = plt.colorbar(cax=position)
        cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

        for t in cbar.ax.get_yticklabels():
            t.set_horizontalalignment('right')
            t.set_x(4.0)

        plt.savefig(out_file + '-chain' + myList[j] + '.png', dpi=200)

        # plt.show()
        plt.close('all')


def interChainCorrelationMaps(ccMatrix, 
                              minColorBarLimit, maxColorBarLimit, 
                              out_file, title, selectedAtoms, saveMatrix):
    """
        Plot inter-chain correlations to different png files,
        if there are at least two chains!

    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    minColorBarLimit: signed int
        Mostly, -1 or 0.
    maxColorBarLimit: unsigned int
        Mostly, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    saveMatrix: bool
        If True, an output file for the correlations will be written
        be written. 

    Returns
    -------
    Nothing
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    myList = list(Counter(selectedAtoms.getChids()).keys())
    # print(myList)
    n = len(ccMatrix)
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

    for k in range(len(myList)):
        for l in range(0, k):
            # Set labels for X axis
            major_numsX = []
            major_labelsX = []
            realTicsListX = np.linspace(selection_reorder[l],
                                        selection_reorder[l + 1], 6, dtype=int)

            major_numsX = np.linspace(selection_reorder[l] - selection_reorder[l],
                                      selection_reorder[l + 1] - selection_reorder[l],
                                      6, dtype=int)

            realTicsListX[-1] = realTicsListX[-1] - 1
            for item in realTicsListX:
                # print(selectedAtoms.getResnums()[item])
                major_labelsX.append(str(selectedAtoms.getResnums()[item]))

            # Set labels for Y axis
            major_numsY = []
            major_labelsY = []
            realTicsListY = np.linspace(selection_reorder[k],
                                        selection_reorder[k + 1], 6, dtype=int)
            major_numsY = np.linspace(selection_reorder[k] - selection_reorder[k],
                                      selection_reorder[k + 1] - selection_reorder[k],
                                      6, dtype=int)

            realTicsListY[-1] = realTicsListY[-1] - 1
            for item in realTicsListY:
                # print(selectedAtoms.getResnums()[item])
                major_labelsY.append(str(selectedAtoms.getResnums()[item]))

            # print(major_numsX)
            # print(major_labelsX)
            ##########################################################################
            # Set plotting parameters
            ##########################################################################
            # Set residue interface definitions
            fig = plt.figure()
            fig.set_size_inches(8.0, 5.5, forward=True)
            plt.rcParams['font.size'] = 16
            ax = fig.add_subplot(1, 1, 1)

            plt.xlabel('Residue indices - Chain ' + myList[l])
            plt.ylabel('Residue indices - Chain ' + myList[k])
            plt.title('Chains ' + myList[l] + '-' + myList[k], y=1.08)

            # plt.rcParams['axes.titlepad'] = 20
            ax.autoscale(False)
            ax.set_aspect('equal')
            # for item in ax.get_xticks():
            #    print(item)

            # ax.set_xticks(major_numsX, major_labelsX)
            plt.xticks(major_numsX, major_labelsX, size=12)
            plt.yticks(major_numsY, major_labelsY, size=12)

            # ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
            # ax.yaxis.set_tick_params(width=2, length=5)

            ax.tick_params(which='major', width=2, length=5)
            ax.tick_params(which='minor', width=1, length=3)

            # print("Min value")
            # print("Row Index of min value")
            # print(np.argmin(ccMatrix_sub, 0))

            # print("Column Index of min value")
            # print(np.argmin(ccMatrix_sub, 1))

            # Set colorbar features here!
            jet = plt.get_cmap('jet')
            djet = cmap_discretize(jet, 8)
            plt.axis([0, n, 0, n])
            plt.axis([0, (selection_reorder[l + 1] - selection_reorder[l]), 0,
                      (selection_reorder[k + 1] - selection_reorder[k])])

            sub_nDCC_matrix = ccMatrix[(selection_reorder[k]): (selection_reorder[k + 1]),
                              (selection_reorder[l]): (selection_reorder[l + 1])]
            plt.imshow(np.matrix(sub_nDCC_matrix), cmap=djet)
            plt.clim(minColorBarLimit, maxColorBarLimit)

            # position=fig.add_axes([0.85, 0.15, 0.03, 0.70])

            divider = make_axes_locatable(plt.gca())
            position = divider.append_axes("right", "5%", pad="25%")
            cbar = plt.colorbar(cax=position)
            cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

            for t in cbar.ax.get_yticklabels():
                t.set_horizontalalignment('right')
                t.set_x(4.0)

            plt.savefig(out_file + '-chains' + myList[l] + '-' + myList[k] + '.png', dpi=200)
            plt.close('all')
            # plt.tight_layout()
            # plt.show()


def distanceDistribution(ccMatrix, out_file, title, selectedAtoms,
                         absoluteValues: bool, writeAllOutput: bool):
    """
        Plot inter-chain correlations vs distances to a png files.

    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    minColorBarLimit: signed int
        Mostly, -1 or 0.
    maxColorBarLimit: unsigned int
        Mostly, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    absoluteValues: bool
        If True, an absolute values of correlations will be consideered. 
    writeAllOutput: bool
        If True, an output file for distances and correlations. 
        This can be useful to see their distribution as well as 
        individual values.  

    Returns
    -------
    Nothing
    """
    # Calculate distance matrix
    dist_matrix = buildDistMatrix(selectedAtoms)

    # Plot the figure
    # print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    x = dist_matrix.flatten()
    y = ccMatrix.flatten()

    # fig, ax = plt.subplots()
    plt.subplots()
    plt.locator_params(axis='y', nbins=4)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Distance ($\AA$)", fontsize=20)
    plt.ylabel(title, fontsize=20)
    if absoluteValues:
        plt.ylim([0.0, 1.0])
        dst_file = out_file + '-absolute-correlation-vs-distance'

    else:
        plt.ylim([-1.0, 1.0])
        plt.axhline(0, color='k', lw=0.5)
        dst_file = out_file + '-correlation-vs-distance'

    plt.plot(x, y, '.', color='k')
    plt.tight_layout()
    # plt.xlim(xmin=0) #xmin will be removed in matplotlib 3.2
    plt.xlim(left=0)
    # plt.show()
    plt.savefig(dst_file + '.png')
    plt.close('all')

    # Write output
    # Writing the output is very important for further analyses such as
    # inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    if writeAllOutput:
        DATA_FILE = open(dst_file + '.dat', 'w')
        for i in range(0, len(ccMatrix)):
            for j in range(i + 1, len(ccMatrix)):
                DATA_FILE.write("{0:d}\t{1:s}\t{2:d}\t{3:s}\t{4:.3f}\t{5:.3f}\n".format(selectedAtoms.getResnums()[i],
                                                                                        selectedAtoms.getChids()[i],
                                                                                        selectedAtoms.getResnums()[j],
                                                                                        selectedAtoms.getChids()[j],
                                                                                        dist_matrix[i][j],
                                                                                        ccMatrix[i][j]))
        DATA_FILE.close()


def projectCorrelationsOntoProteinVMD(pdb_file, ccMatrix, vmd_out_file,
                                      selectedAtoms, valueFilter,
                                      absoluteValues: bool,
                                      writeAllOutput: bool):
    """
        This function writes tcl files that contains the correlations between
    residues i and j. It produces three output files:
    1-A general file that contains all correlation.
    2-(If there are at least two chains) Files that contain interchain
      correlations.
    3-(If there are at least two chains) Files that contain intrachain
      correlations of individual chains.
    The output files can be visualized with VMD (Visual Molecular
    dynamics) program as follows.
    i) Load your pdb file, whether via command line or graphical interface.
    ii) Go to Extemsions -> Tk Console and then
    iii) source vmd-output-general.tcl
    It can take some to load the general script.

    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    vmd_out_file: string
        prefix for the output tcl files.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    valueFilter: float
        Correlation values smaller than this threshold will not be written
        for visualization. For example, 0.3 is a good threshold for normalized
        dynamical cross-correlation data. 
    absoluteValues: bool
        If True, an absolute values of correlations will be consideered. 
    writeAllOutput: bool
        If True, an output file for distances and correlations. 
        This can be useful to see their distribution as well as 
        individual values.  

    Returns
    -------
    Nothing
    """
    # Calculate distance matrix
    dist_matrix = buildDistMatrix(selectedAtoms)

    # Plot the figure
    # print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    x = dist_matrix.flatten()
    y = ccMatrix.flatten()

    # print(len(y))

    distanceFilter = 0.5

    # Write output in VMD format
    # Writing the output is very important for further analyses such as
    # inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    #
    draw_string = "draw cylinder " \
                  " [lindex [[atomselect top \"chain {0:s} and resid {1:d} and name CA\"] get {{x y z}}] 0]" \
                  " [lindex [[atomselect top \"chain {2:s} and resid {3:d} and name CA\"] get {{x y z}}] 0]" \
                  " radius {4:.3f}\n"
    vdw_representation_string = "mol representation VDW 0.750000 50.000000\n" + \
                                "mol material Glossy\n" + \
                                "#mol color ColorID 7\n" + \
                                "mol selection \"chain {0:s} and resid {1:d} and name CA\"\n" + \
                                "mol addrep 0\n"
    DATA_FILE = open(vmd_out_file + '-general.tcl', 'w')
    DATA_FILE.write("mol new " + pdb_file + " \n")
    # DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
    DATA_FILE.write("mol modstyle 0 0 Tube\n")
    DATA_FILE.write("mol modcolor 0 0 Chain\n")
    DATA_FILE.write("mol modmaterial 0 0 MetallicPastel\n")
    for i in range(0, len(ccMatrix)):
        for j in range(i + 1, len(ccMatrix)):
            if (np.absolute(ccMatrix[i][j]) > valueFilter):
                DATA_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[i],
                                                                 selectedAtoms.getResnums()[i]))
                DATA_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[j],
                                                                 selectedAtoms.getResnums()[j]))
                DATA_FILE.write(draw_string.format(selectedAtoms.getChids()[i],
                                                   selectedAtoms.getResnums()[i],
                                                   selectedAtoms.getChids()[j],
                                                   selectedAtoms.getResnums()[j],
                                                   # The radius of the connecting cylinder is proportional
                                                   # to the correlation value.
                                                   # However, it is necessary to multiply the radius
                                                   # with 0.5 to make it look better.
                                                   ccMatrix[i][j] * 0.5))
    DATA_FILE.close()

    chains = Counter(selectedAtoms.getChids()).keys()

    plotChains = True
    if (len(chains) > 1) & plotChains:
        # Inter-chain
        for chainI in chains:
            for chainJ in chains:
                if chainI != chainJ:
                    DATA_FILE = open(vmd_out_file + '-interchain-chains' + chainI + '-' + chainJ + '.tcl', 'w')
                    DATA_FILE.write("mol new " + pdb_file + " \n")
                    #DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
                    DATA_FILE.write("mol modstyle 0 0 Tube\n")
                    DATA_FILE.write("mol modcolor 0 0 Chain\n")
                    DATA_FILE.write("mol modmaterial 0 0 MetallicPastel\n")
                    for i in range(0, len(ccMatrix)):
                        for j in range(i + 1, len(ccMatrix)):
                            if np.absolute(ccMatrix[i][j]) > valueFilter:
                                if (selectedAtoms.getChids()[i] == chainI) and (selectedAtoms.getChids()[j] == chainJ):
                                    DATA_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[i],
                                                                                     selectedAtoms.getResnums()[i]))
                                    DATA_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[j],
                                                                                     selectedAtoms.getResnums()[j]))

                                    DATA_FILE.write(draw_string.format(selectedAtoms.getChids()[i],
                                                                       selectedAtoms.getResnums()[i],
                                                                       selectedAtoms.getChids()[j],
                                                                       selectedAtoms.getResnums()[j], \
                                                                       # The radius of the connecting cylinder is
                                                                       # proportional to the correlation value.
                                                                       # However, it is necessary to multiply
                                                                       # the radius with 0.5 to make it look better.
                                                                       ccMatrix[i][j] * 0.5))
                    DATA_FILE.close()

        # Intra-chain
        for chain in chains:
            DATA_FILE = open(vmd_out_file + '-intrachain-chain' + chain + '.tcl', 'w')
            DATA_FILE.write("mol new " + pdb_file + " \n")
            #DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
            DATA_FILE.write("mol modstyle 0 0 Tube\n")
            DATA_FILE.write("mol modcolor 0 0 Chain\n")
            DATA_FILE.write("mol modmaterial 0 0 MetallicPastel\n")
            for i in range(0, len(ccMatrix)):
                for j in range(i + 1, len(ccMatrix)):
                    if np.absolute(ccMatrix[i][j]) > valueFilter:
                        if (selectedAtoms.getChids()[i] == chain) and (selectedAtoms.getChids()[j] == chain):
                            DATA_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[i],
                                                                             selectedAtoms.getResnums()[i]))
                            DATA_FILE.write(vdw_representation_string.format(selectedAtoms.getChids()[j],
                                                                             selectedAtoms.getResnums()[j]))
                            DATA_FILE.write(draw_string.format(selectedAtoms.getChids()[i],
                                                               selectedAtoms.getResnums()[i],
                                                               selectedAtoms.getChids()[j],
                                                               selectedAtoms.getResnums()[j],
                                                               # The radius of the connecting cylinder is proportional
                                                               # to the correlation value.
                                                               # However, it is necessary to multiply the radius
                                                               # with 0.5 to make it look better.
                                                               ccMatrix[i][j] * 0.5))
            DATA_FILE.close()


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
    numOfLabels=9
    generatePNG(diffMap, minColorBarLimit, maxColorBarLimit, numOfLabels,
                            out_file+"-overall-difference.png", title, selectedAtoms)

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

def triangulateMaps(ccMatrix1, ccMatrix2,
                                minColorBarLimit, maxColorBarLimit, 
                                out_file, title, selectedAtoms):
    """
        Given two correlation maps, it puts them into upper and lower triangles.
        Sizes of ccMatrix1 and ccMatrix2 are identical. Only one atom set is
        sufficient to identify the residue IDs.
        
    Parameters
    ----------
    ccMatrix1: A numpy square matrix of floats
        The first correlation matrix.
    ccMatrix2: A numpy square matrix of floats
        The second correlation matrix.
    minColorBarLimit: signed int
        If nDCC maps -1. If absndcc or lmi, 0.  
    maxColorBarLimit: unsigned int
        If nDCC maps 1. If absndcc or lmi, 1.
    out_file: string
        prefix for the output png files.
        This prefix will get _overall.png extension.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.

    Returns
    -------
    ccMatrixCombined: A numpy square matrix of floats
        This is a matrix that contain ccMatrix1 in the upper triangle
        and the ccMatrix2 in the lower triangle. 
    """

    n = len(ccMatrix1)
    if (len(ccMatrix2) != n):
        print("Two matrices have to have same dimensions!")
        sys.exit(-1)

    ccMatrixCombined = np.zeros((n, n), np.double)
    ccMatrixCombined = np.triu(ccMatrix2, k=0) + np.tril(ccMatrix1, k=0)                                                                                                                                             
    np.fill_diagonal(ccMatrixCombined, 1.0)
    
    numOfLabels=5
    generatePNG(ccMatrixCombined, minColorBarLimit, maxColorBarLimit,
                numOfLabels, out_file+"-triangulated.png", title, selectedAtoms)
    return ccMatrixCombined

def generatePNG(ccMatrix, minColorBarLimit, maxColorBarLimit,
                numOfLabels, out_file, title, selectedAtoms):
    """
        Generates a heatmap in PNG format from the ccMatrix. 
        
    Parameters
    ----------
    ccMatrix: A numpy square matrix of floats
        The first correlation matrix.
    minColorBarLimit: signed int
        If nDCC maps -1. If absndcc or lmi, 0.  
    maxColorBarLimit: unsigned int
        If nDCC maps 1. If absndcc or lmi, 1.
    numOfLabels: int
        Number of labels on colorbar. 
    out_file: string
        prefix for the output png files.
    title: string
        Title of the figure.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.

    Returns
    -------
    Nothing 
    """
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
    n = len(ccMatrix)
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

    plt.imshow(np.matrix(ccMatrix), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position = fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar = plt.colorbar(cax=position, format='%.2f')

    #cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])
    
    cbar.set_ticks(np.linspace(minColorBarLimit, maxColorBarLimit, num=numOfLabels))
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
    plt.savefig(out_file, bbox_inches='tight', dpi=200)
    # plt.show()


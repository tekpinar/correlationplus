##############################################################################
# correlationPlus - Python module to plot dynamical correlations maps         #
#                   for proteins.                                             #
# Authors: Mustafa Tekpinar                                                   #
# Copyright Mustafa Tekpinar 2017-2018                                        #
# Copyright CNRS-UMR3528, 2019                                                #
# Copyright Institut Pasteur Paris, 2020                                       #
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

from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
from prody import buildDistMatrix


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
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
    data_file = open(inp_file, 'r')

    allLines = data_file.readlines()
    data_list = []
    for line in allLines:
        words = line.split()
        for i in words:
                data_list.append(i)

    data_list = data_list[4:-1]
    n = int(np.sqrt(len(data_list)))
    data_array = np.array(data_list, dtype=float)

    cc = np.reshape(data_array, (n, n))
    #print (cc.dtype)

    data_file.close()
    #Fill diagonal elements with zeros
    #np.fill_diagonal(cc, 0.0)

    #Find maximum for the remaining matrix
    #maximum=cc.max()
    #cc=cc/cc.max()
    #print (maximum)
    #np.fill_diagonal(cc, 1.0)

    if writeAllOutput:
        np.savetxt(inp_file[:-4] + "_modif.dat", cc, fmt='%.6f')
    return cc


def filterCorrelationMapByDistance(ccMatrix, out_file, title,
                                   selectedAtoms, distanceValue,
                                   absoluteValues: bool, writeAllOutput: bool):
    """
    If residues are closer to each other than a certain distance
    (distanceValue), make these correlations zero. This filter can be useful
    to get rid of short distance correlation for visualization purposes.
    This function returns a filtered ccMatrix.
    """
    print("@> Filtering correlations lower than " + str(distanceValue) + " Angstrom")
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
    minColorBarLimit: unsigned int
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
    Plot intra-chain correlations if there are at least two chains!
    """

    myList = list(Counter(selectedAtoms.getChids()).keys())
    # print(myList)

    # selection_tick_labels=[]
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


def interChainCorrelationMaps(ccMatrix, minColorBarLimit,
                              maxColorBarLimit, out_file, title, selectedAtoms, saveMatrix):
    """
    Plot inter-chain correlations if there are at least two chains!
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    myList = list(Counter(selectedAtoms.getChids()).keys())
    # print(myList)
    n = len(ccMatrix)
    # selection_tick_labels=[]
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
    # Calculate distance matrix
    dist_matrix = buildDistMatrix(selectedAtoms)

    # Plot the figure
    # print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    x = dist_matrix.flatten()
    y = ccMatrix.flatten()

    # print(len(y))

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


def projectCorrelationsOntoProteinVMD(ccMatrix, vmd_out_file,
                                      selectedAtoms, valueFilter,
                                      absoluteValues: bool,
                                      writeAllOutput: bool):
    """
    This function writes tcl files that contains the correlations between
    residues i and j. It produces three output files:
    1-A general file that contains all correlation.
    2-(If there are at least two chains) A file that contains interchain
      correlations.
    3-(If there are at least two chains) A file that contains intrachain
      correlations.
    The output files can be visualized with VMD (Visual Molecular
    dynamics) program as follows.
    i) Load your pdb file, whether via command line or graphical interface.
    ii) Go to Extemsions -> Tk Console and then
    iii) source vmd-output-general.tcl
    It can take some to load the general script.
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
    # valueFilter = 0.5
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

    DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
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
                    DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
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
            DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
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

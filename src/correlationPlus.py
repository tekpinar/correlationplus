#!/usr/bin/env python
"""
Program Name: Cross-correlation Plotting Program (I'll find a fancy name later!)
Author      : Mustafa TEKPINAR
Email       : tekpinar@buffalo.edu
Copyright   : Mustafa Tekpinar - 2020
License     : MIT License

Purpose     : This is a small program to automatize plotting of normalized 
dynamical cross-correlations obtained from molecular dynamics simulations or 
elastic network models. This script can be useful if you have multiple 
chains in a structure and you want to see intra-chain and inter-chain 
correlations more clearly. I just didn't like the way current programs are 
doing it and I wrote something for myself. I hope it may help the others also!
"""

import matplotlib.pyplot as plt
import matplotlib.collections as collections
import numpy as np
import getopt
import sys
import matplotlib

#from prody import *
from prody import parsePDB
from prody import buildDistMatrix
from collections import Counter, OrderedDict


import networkx as nx
from math import fabs
from math import log

def usage_correlationMaps():
    """
    Show how to use this program!
    """
    print("\nExample usage:\n")
    print("python correlationPlus.py -i 4z90-cross-correlations.txt -p 4z90.pdb \n")
    print("Arguments: -i: A file containing normalized dynamical cross correlations in matrix format. (Mandatory)")
    print("           -p: PDB file of the protein. (Mandatory)")
    print("           -s: It can be ndcc, lmi or absndcc (absolute values of ndcc). Default value is ndcc (Optional)")
    print("           -o: This will be your output file. Output figures are in png format. (Optional)\n\n")

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
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

##################### HANDLE ARGUMENTS ##########################################                                                                           
def handle_arguments_correlationMaps():
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:s:p:",["help", "inp=", "out=", "sel=", "pdb="])
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:s:p:",["help", "inp=", "out=", "sel=", "pdb="])
    except getopt.GetoptError:
        usage_correlationMaps()
    for opt, arg in opts:
        # if opt == '-h':
        #     usage_correlationMaps()
        #     sys.exit(-1)
        if opt in ('-h', "--help"):
            usage_correlationMaps()
            sys.exit(-1)
        elif opt in ("-i", "--inp"):
            inp_file = arg
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-s", "--sel"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file = arg
        else:
            assert False, usage_correlationMaps()

    #Input data matrix and PDB file are mandatory!
    if inp_file==None or pdb_file==None:
        usage_correlationMaps()
        sys.exit(-1)

    #Assign a default name if the user forgets the output file prefix.
    if (out_file == None):
        out_file = "correlation"

    #The user may prefer not to submit a title for the output.
    if (sel_type == None):
        sel_type = "ndcc"

    return (inp_file, out_file, sel_type, pdb_file)

def overallCorrelationMap(ccMatrix, minColorBarLimit, maxColorBarLimit, out_file, title, selectedAtoms):
    """
    Plots nDCC maps for the whole structure
    """

    #selectedAtoms = parsePDB(pdb_file, subset='ca')
    #chainLengths = Counter(selectedAtoms.getChids()).values()

    n=len(ccMatrix)
    ##########################################################################
    #Set residue interface definitions
    fig=plt.figure()
    fig.set_size_inches(8.0, 5.5, forward=True)
    plt.rcParams['font.size'] = 16
    ax=fig.add_subplot(1,1,1)

    plt.xlabel('Residue indices')
    plt.ylabel('Residue indices')
    plt.title(title, y=1.08)

    #print(selectedAtoms.getChids())

    

    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtoms.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    myList = list(Counter(selectedAtoms.getChids()).keys())

    major_nums=[]
    major_labels=[]

    if(len(myList)==1):
        realTicsList = np.linspace(0, len(selectedAtoms.getResnums())-1, 6, dtype=int)
        major_nums = realTicsList

        for item in (realTicsList):
            #print(selectedAtoms.getResnums()[item])
            major_labels.append(str(selectedAtoms.getResnums()[item]))
    elif(len(myList)>1):
        for i in Counter(selectedAtoms.getChids()).values():
            tempVal = tempVal + i
            selection_reorder.append(tempVal)
            selection_tick_labels.append(str(selectedAtoms.getResnums()[tempVal-1]))
        major_nums.extend(selection_reorder)
        major_labels.extend(selection_tick_labels)
    else:
        print("Warning: Unknown chain ID!")
    #print(selection_reorder)

    ##########################################################################
    #Set plotting parameters

    #plt.rcParams['axes.titlepad'] = 20
    ax.autoscale(False)
    ax.set_aspect('equal')

    #ax.set_xticks(major_nums, major_labels, rotation=45, minor=False)
    plt.xticks(major_nums, major_labels, size=12, rotation=45)
    plt.yticks(major_nums, major_labels, size=12)

    #ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
    #ax.yaxis.set_tick_params(width=2, length=5)

    plt.axis([0, n, 0, n])
    ax.tick_params(which='major', width=2, length=5)
    ax.tick_params(which='minor', width=1, length=3)

    #print("Min value")
    #print("Row Index of min value")
    #print(np.argmin(ccMatrix_sub, 0))

    #print("Column Index of min value")
    #print(np.argmin(ccMatrix_sub, 1))

    #Set colorbar features here!  
    jet=plt.get_cmap('jet') 
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(ccMatrix), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position=fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar=plt.colorbar(cax=position)

    cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')   
        t.set_x(4.0)

    if(len(myList)>1):
        #Add chain borders to the plot
        for i in range(len(selection_reorder)-1):
            beginningPoint = selection_reorder[i]/selection_reorder[-1]
            endingPoint = selection_reorder[i+1]/selection_reorder[-1]
            middlePoint = (float(beginningPoint)+float(endingPoint))/2.0
            if(i%2==0):
                #x axis
                ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                            xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black'))

                #y axis
                ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                            xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black'))  

            elif(i%2==1):
                #x axis
                ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                            xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray'))
                
                #y axis
                ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                            xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray')) 

            ax.annotate(myList[i], xy=(0, 1.04), xycoords='axes fraction', xytext=(middlePoint-0.015, 1.04), size=14, color='black')
            ax.annotate(myList[i], xy=(1.05, 0), xycoords='axes fraction', xytext=(1.05, middlePoint-0.015), rotation=90, size=14, color='black')
            #print(middlePoint)

    #plt.tight_layout()       
    plt.savefig(out_file+'-overall.png', bbox_inches='tight', dpi=200)
    #plt.show()

def intraChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit, out_file, title, selectedAtoms, saveMatrix):
    """
    Plot intra-chain correlations if there are at least two chains!
    """

    myList = list(Counter(selectedAtoms.getChids()).keys())
    #print(myList)

    #selection_tick_labels=[]
    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtoms.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    for i in Counter(selectedAtoms.getChids()).values():
        tempVal = tempVal + i
        selection_reorder.append(tempVal)
        selection_tick_labels.append(str(selectedAtoms.getResnums()[tempVal-1]))

    #print(selection_reorder)

    for j in range(len(myList)):
        #Set labels
        major_nums=[]
        major_labels=[]
        realTicsList = np.linspace(selection_reorder[j], selection_reorder[j+1], 6, dtype=int)
        ticsList = np.linspace(selection_reorder[j]-selection_reorder[j], selection_reorder[j+1]-selection_reorder[j], 6, dtype=int)
        #selection_reorder[j+1]-
        major_nums = ticsList
        #major_nums = major_nums[0:-1]
        #print(realTicsList)
        realTicsList[-1] = realTicsList[-1]-1
        for item in realTicsList:
            #print(selectedAtoms.getResnums()[item])
            major_labels.append(str(selectedAtoms.getResnums()[item]))
        
        #print(major_nums)
        #print(major_labels)
        ##########################################################################
        #Set plotting parameters
        ##########################################################################
        #Set residue interface definitions
        fig=plt.figure()
        fig.set_size_inches(8.0, 5.5, forward=True)
        plt.rcParams['font.size'] = 16
        ax=fig.add_subplot(1,1,1)

        plt.xlabel('Residue indices')
        plt.ylabel('Residue indices')
        plt.title('Chain '+myList[j], y=1.08)

        #plt.rcParams['axes.titlepad'] = 20
        ax.autoscale(False)
        ax.set_aspect('equal')
        #for item in ax.get_xticks():
        #    print(item)

        #ax.set_xticks(major_nums, major_labels)
        plt.xticks(major_nums, major_labels, size=12)
        plt.yticks(major_nums, major_labels, size=12)

        #ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
        #ax.yaxis.set_tick_params(width=2, length=5)

        ax.tick_params(which='major', width=2, length=5)
        ax.tick_params(which='minor', width=1, length=3)

        #print("Min value")
        #print("Row Index of min value")
        #print(np.argmin(ccMatrix_sub, 0))

        #print("Column Index of min value")
        #print(np.argmin(ccMatrix_sub, 1))

        #Set colorbar features here!
        jet=plt.get_cmap('jet')
        djet = cmap_discretize(jet, 8)
        plt.axis([0, (selection_reorder[j+1]-selection_reorder[j]), 0, (selection_reorder[j+1]-selection_reorder[j])])
        sub_nDCC_matrix = ccMatrix[(selection_reorder[j]) : (selection_reorder[j+1]), (selection_reorder[j]) : (selection_reorder[j+1])]
        plt.imshow(np.matrix(sub_nDCC_matrix), cmap=djet)
        plt.clim(minColorBarLimit, maxColorBarLimit)

        position=fig.add_axes([0.85, 0.15, 0.03, 0.70])
        
        cbar=plt.colorbar(cax=position)
        cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

        for t in cbar.ax.get_yticklabels():
            t.set_horizontalalignment('right')   
            t.set_x(4.0)
        
        plt.savefig(out_file+'-chain'+myList[j]+'.png', dpi=200)

        #plt.show()
        plt.close('all')

def interChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit, out_file, title, selectedAtoms, saveMatrix):
    """
    Plot inter-chain correlations if there are at least two chains!
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    myList = list(Counter(selectedAtoms.getChids()).keys())
    #print(myList)
    n=len(ccMatrix)
    #selection_tick_labels=[]
    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtoms.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    for i in Counter(selectedAtoms.getChids()).values():
        tempVal = tempVal + i
        selection_reorder.append(tempVal)
        selection_tick_labels.append(str(selectedAtoms.getResnums()[tempVal-1]))

    #print(selection_reorder)

    for k in range(len(myList)):
        for l in range(0, k):
            #Set labels for X axis
            major_numsX=[]
            major_labelsX=[]
            realTicsListX = np.linspace(selection_reorder[l], \
                                        selection_reorder[l+1], 6, dtype=int)

            major_numsX = np.linspace(selection_reorder[l]-selection_reorder[l], \
                                    selection_reorder[l+1]-selection_reorder[l],\
                                                                6, dtype=int)

            realTicsListX[-1] = realTicsListX[-1]-1
            for item in realTicsListX:
                #print(selectedAtoms.getResnums()[item])
                major_labelsX.append(str(selectedAtoms.getResnums()[item]))

            #Set labels for Y axis
            major_numsY=[]
            major_labelsY=[]
            realTicsListY = np.linspace(selection_reorder[k], \
                                        selection_reorder[k+1], 6, dtype=int)
            major_numsY = np.linspace(selection_reorder[k]-selection_reorder[k], \
                                    selection_reorder[k+1]-selection_reorder[k], 
                                                                6, dtype=int)
            

            realTicsListY[-1] = realTicsListY[-1]-1
            for item in realTicsListY:
                #print(selectedAtoms.getResnums()[item])
                major_labelsY.append(str(selectedAtoms.getResnums()[item]))
            
            #print(major_numsX)
            #print(major_labelsX)
            ##########################################################################
            #Set plotting parameters
            ##########################################################################
            #Set residue interface definitions
            fig=plt.figure()
            fig.set_size_inches(8.0, 5.5, forward=True)
            plt.rcParams['font.size'] = 16
            ax=fig.add_subplot(1,1,1)

            plt.xlabel('Residue indices - Chain '+myList[l])
            plt.ylabel('Residue indices - Chain '+myList[k])
            plt.title('Chains '+myList[l]+'-'+myList[k], y=1.08)

            #plt.rcParams['axes.titlepad'] = 20
            ax.autoscale(False)
            ax.set_aspect('equal')
            #for item in ax.get_xticks():
            #    print(item)

            #ax.set_xticks(major_numsX, major_labelsX)
            plt.xticks(major_numsX, major_labelsX, size=12)
            plt.yticks(major_numsY, major_labelsY, size=12)

            #ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
            #ax.yaxis.set_tick_params(width=2, length=5)

            ax.tick_params(which='major', width=2, length=5)
            ax.tick_params(which='minor', width=1, length=3)

            #print("Min value")
            #print("Row Index of min value")
            #print(np.argmin(ccMatrix_sub, 0))

            #print("Column Index of min value")
            #print(np.argmin(ccMatrix_sub, 1))

            #Set colorbar features here!
            jet=plt.get_cmap('jet')
            djet = cmap_discretize(jet, 8)
            plt.axis([0, n, 0, n])
            plt.axis([0, (selection_reorder[l+1]-selection_reorder[l]), 0, \
                         (selection_reorder[k+1]-selection_reorder[k])])

            sub_nDCC_matrix = ccMatrix[(selection_reorder[k]) : (selection_reorder[k+1]), \
                                       (selection_reorder[l]) : (selection_reorder[l+1])]
            plt.imshow(np.matrix(sub_nDCC_matrix), cmap=djet)
            plt.clim(minColorBarLimit, maxColorBarLimit)

            #position=fig.add_axes([0.85, 0.15, 0.03, 0.70])
            
            divider = make_axes_locatable(plt.gca())
            position = divider.append_axes("right", "5%", pad="25%")
            cbar=plt.colorbar(cax=position)
            cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

            for t in cbar.ax.get_yticklabels():
                t.set_horizontalalignment('right')   
                t.set_x(4.0)
            
            plt.savefig(out_file+'-chains'+myList[l]+'-'+myList[k]+'.png', dpi=200)
            plt.close('all')
            #plt.tight_layout()
            #plt.show()

def distanceDistribution(ccMatrix, out_file, title, selectedAtoms, \
    absoluteValues: bool, writeAllOutput: bool):
    #Calculate distance matrix
    dist_matrix=buildDistMatrix(selectedAtoms)
    
    #Plot the figure
    #print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    x=dist_matrix.flatten()
    y=ccMatrix.flatten()

    #print(len(y))

    #fig, ax = plt.subplots()
    plt.subplots()
    plt.locator_params(axis='y', nbins=4)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Distance ($\AA$)", fontsize=20)
    plt.ylabel(title, fontsize=20)
    if(absoluteValues):
        plt.ylim([0.0, 1.0])
        dst_file = out_file+'-absolute-correlation-vs-distance'

    else:
        plt.ylim([-1.0, 1.0])
        plt.axhline(0,color='k',lw=0.5)
        dst_file = out_file+'-correlation-vs-distance'

    plt.plot(x,y, '.', color='k')
    plt.tight_layout()
    plt.xlim(xmin=0)
    #plt.show()
    plt.savefig(dst_file+'.png')
    plt.close('all')
    
    #Write output
    #Writing the output is very important for further analyses such as 
    #inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    if(writeAllOutput):
        DATA_FILE = open(dst_file+'.dat', 'w')
        for i in range(0, len(ccMatrix)):
            for j in range(i+1, len(ccMatrix)):
                DATA_FILE.write("{0:d}\t{1:s}\t{2:d}\t{3:s}\t{4:.3f}\t{5:.3f}\n".\
                                format(selectedAtoms.getResnums()[i],\
                                        selectedAtoms.getChids()[i],\
                                        selectedAtoms.getResnums()[j],\
                                        selectedAtoms.getChids()[j],\
                                        dist_matrix[i][j],\
                                        ccMatrix[i][j]))
        DATA_FILE.close()
    
def convertLMIdata2Matrix(inp_file, writeAllOutput: bool):
    data_file=open(inp_file, 'r')

    allLines=data_file.readlines()
    data_list=[]
    for line in allLines:
        words=line.split()
        for i in words:
                data_list.append(i)

    data_list = data_list[4:-1]
    n=int(np.sqrt(len(data_list)))
    data_array=np.array(data_list, dtype=float)

    cc=np.reshape(data_array, (n, n))
    #print (cc.dtype)

    data_file.close()
    #Fill diagonal elements with zeros
    #np.fill_diagonal(cc, 0.0)

    #Find maximum for the remaining matrix
    #maximum=cc.max()
    #cc=cc/cc.max()
    #print (maximum)
    #np.fill_diagonal(cc, 1.0)

    if(writeAllOutput):
        np.savetxt(inp_file[:-4]+"_modif.dat", cc, fmt='%.6f')
    return cc

# def readForceConstantsMatrix(inp_file, writeAllOutput: bool):
# """
#     The force constant matrix is in the following format!
#     i j k_ij
# """
#     data_file=open(inp_file, 'r')

#     allLines=data_file.readlines()
#     data_list=[]
#     for line in allLines:
#         words=line.split()
#         for i in words:
#                 data_list.append(i)

#     #data_list = data_list[4:-1]
#     n=int(np.sqrt(len(data_list)))
#     data_array=np.array(data_list, dtype=float)

#     cc=np.reshape(data_array, (n, n))
#     #print (cc.dtype)

#     data_file.close()
#     #Fill diagonal elements with zeros
#     #np.fill_diagonal(cc, 0.0)

#     #Find maximum for the remaining matrix
#     #maximum=cc.max()
#     #cc=cc/cc.max()
#     #print (maximum)
#     #np.fill_diagonal(cc, 1.0)

#     if(writeAllOutput):
#         np.savetxt(inp_file[:-4]+"_modif.dat", cc, fmt='%.6f')
#     return cc

def overallUniformDifferenceMap(ccMatrix1, ccMatrix2, minColorBarLimit, maxColorBarLimit, out_file, title, selectedAtoms):
    """
    Plots the difference map between correlation maps for the entire structure. 
    Sizes of ccMatrix1 and ccMatrix2 are identical. Only one atom set is 
    sufficient to plot the difference map. 
    """

    #selectedAtoms = parsePDB(pdb_file, subset='ca')
    #chainLengths = Counter(selectedAtoms.getChids()).values()
    diffMap = np.subtract(ccMatrix1, ccMatrix2)

    n=len(ccMatrix1)
    ##########################################################################
    #Set residue interface definitions
    fig=plt.figure()
    fig.set_size_inches(8.0, 5.5, forward=True)
    plt.rcParams['font.size'] = 16
    ax=fig.add_subplot(1,1,1)

    plt.xlabel('Residue indices')
    plt.ylabel('Residue indices')
    plt.title(title, y=1.08)
    plt.grid(color='w', linestyle='--', linewidth=1)

    #print(selectedAtoms.getChids())

    myList = list(Counter(selectedAtoms.getChids()).keys())

    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtoms.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    for i in Counter(selectedAtoms.getChids()).values():
        tempVal = tempVal + i
        selection_reorder.append(tempVal)
        selection_tick_labels.append(str(selectedAtoms.getResnums()[tempVal-1]))

    #print(selection_reorder)

    major_nums=[]
    major_labels=[]
    major_nums.extend(selection_reorder)
    major_labels.extend(selection_tick_labels)

    ##########################################################################
    #Set plotting parameters

    #plt.rcParams['axes.titlepad'] = 20
    ax.autoscale(False)
    ax.set_aspect('equal')

    #ax.set_xticks(major_nums, major_labels, rotation=45, minor=False)
    plt.xticks(major_nums, major_labels, size=12, rotation=45)
    plt.yticks(major_nums, major_labels, size=12)

    #ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
    #ax.yaxis.set_tick_params(width=2, length=5)

    plt.axis([0, n, 0, n])
    ax.tick_params(which='major', width=2, length=5)
    ax.tick_params(which='minor', width=1, length=3)

    #print("Min value")
    #print("Row Index of min value")
    #print(np.argmin(ccMatrix1_sub, 0))

    #print("Column Index of min value")
    #print(np.argmin(ccMatrix1_sub, 1))

    #Set colorbar features here!  
    jet=plt.get_cmap('jet') 
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(diffMap), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position=fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar=plt.colorbar(cax=position)

    cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')   
        t.set_x(4.0)

    #Add chain borders to the plot
    for i in range(len(selection_reorder)-1):
        beginningPoint = selection_reorder[i]/selection_reorder[-1]
        endingPoint = selection_reorder[i+1]/selection_reorder[-1]
        middlePoint = (float(beginningPoint)+float(endingPoint))/2.0
        if(i%2==0):
            #x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black'))

            #y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black'))  

        elif(i%2==1):
            #x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray'))
            
            #y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray')) 

        ax.annotate(myList[i], xy=(0, 1.04), xycoords='axes fraction', xytext=(middlePoint-0.015, 1.04), size=14, color='black')
        ax.annotate(myList[i], xy=(1.05, 0), xycoords='axes fraction', xytext=(1.05, middlePoint-0.015), rotation=90, size=14, color='black')
        #print(middlePoint)

    #plt.tight_layout()       
    plt.savefig(out_file+'-overall-difference.png', bbox_inches='tight', dpi=200)
    #plt.show()
def findCommonCorePDB(selectedAtomSet1, selectedAtomSet2):
    """
    This function assumes that two structures are obtained 
    exactly from the same species and there is not any mutation in any one
    of them. This can happen when two conformations of a protein are obtained
    with different number of missing atoms. 
    Under these conditions, we are trying to match two structures and find
    the common core of two structures that have different number of CA atoms.
    We will use this information to subtract two correlation maps!
    Here, the output will be a dictonary of indices matching in two structures.
    """
    print("@> Calculating a common core for two structures.")
    lengthSet1 = len(selectedAtomSet1)
    lengthSet2 = len(selectedAtomSet2)
    
    #print(lengthSet1)
    #print(lengthSet2)

    #I made it an ordered dictionary so that the it will not shuffle 
    #the items in the dictionary!
    commonCoreDictionary = OrderedDict()
    # if (lengthSet1 > lengthSet2):
    for i in range (0, lengthSet1):
        for j in range (0, lengthSet2):
            if((selectedAtomSet1.getResnums()[i] == \
                    selectedAtomSet2.getResnums()[j]) and 
                (selectedAtomSet1.getResnames()[i] == \
                    selectedAtomSet2.getResnames()[j]) and \
                (selectedAtomSet1.getChids()[i] == \
                    selectedAtomSet2.getChids()[j])):
                commonCoreDictionary[i] = j
    #print(commonCoreDictionary)
    return commonCoreDictionary

def overallNonUniformDifferenceMap(ccMatrix1, ccMatrix2, minColorBarLimit, \
                                    maxColorBarLimit, out_file, title, \
                                    selectedAtomSet1, selectedAtomSet2):
    """
    Plots the difference map between correlation maps for the entire structure. 
    Sizes of ccMatrix1 and ccMatrix2 are not identical. A mapping for matching 
    residues is performed before difference map plotting. 
    """

    #selectedAtoms = parsePDB(pdb_file, subset='ca')
    #chainLengths = Counter(selectedAtoms.getChids()).values()

    commonCoreDictionary = findCommonCorePDB(selectedAtomSet1, selectedAtomSet2)
    diffMap = np.zeros((len(commonCoreDictionary), len(commonCoreDictionary)), dtype=float)

    items = list(commonCoreDictionary.items())

    for i in range (0, len(commonCoreDictionary)):
        key1, value1 = items[i]
        for j in range (i, len(commonCoreDictionary)):
            key2, value2 = items[j]
            diffMap[i][j] = ccMatrix1[key1][key2]-ccMatrix2[value1][value2]
            diffMap[j][i] = diffMap[i][j]  
    #diffMap = np.subtract(ccMatrix1, ccMatrix2)

    n=len(commonCoreDictionary)
    ##########################################################################
    #Set residue interface definitions
    fig=plt.figure()
    fig.set_size_inches(8.0, 5.5, forward=True)
    plt.rcParams['font.size'] = 16
    ax=fig.add_subplot(1,1,1)

    plt.xlabel('Residue indices')
    plt.ylabel('Residue indices')
    plt.title(title, y=1.08)
    plt.grid(color='w', linestyle='--', linewidth=1)

    ##########################################################################
    #Set xtics, ytic, xlabels, y labels etc. 
    #This part may have to be rewritten! 
    myList = list(Counter(selectedAtomSet1.getChids()).keys())

    selection_reorder = []
    selection_tick_labels = []
    selection_tick_labels.append(str(selectedAtomSet1.getResnums()[0]))
    selection_reorder.append(0)
    tempVal = 0
    for i in Counter(selectedAtomSet1.getChids()).values():
        tempVal = tempVal + i
        selection_reorder.append(tempVal)
        selection_tick_labels.append(str(selectedAtomSet1.getResnums()[tempVal-1]))

    #print(selection_reorder)

    major_nums=[]
    major_labels=[]
    major_nums.extend(selection_reorder)
    major_labels.extend(selection_tick_labels)

    ##########################################################################
    #Set plotting parameters

    #plt.rcParams['axes.titlepad'] = 20
    ax.autoscale(False)
    ax.set_aspect('equal')

    #ax.set_xticks(major_nums, major_labels, rotation=45, minor=False)
    plt.xticks(major_nums, major_labels, size=12, rotation=45)
    plt.yticks(major_nums, major_labels, size=12)

    #ax.xaxis.set_tick_params(width=2, length=5, labelsize=12, minor=False)
    #ax.yaxis.set_tick_params(width=2, length=5)

    plt.axis([0, n, 0, n])
    ax.tick_params(which='major', width=2, length=5)
    ax.tick_params(which='minor', width=1, length=3)

    #print("Min value")
    #print("Row Index of min value")
    #print(np.argmin(ccMatrix1_sub, 0))

    #print("Column Index of min value")
    #print(np.argmin(ccMatrix1_sub, 1))

    #Set colorbar features here!  
    jet=plt.get_cmap('jet') 
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(diffMap), cmap=djet)
    plt.clim(minColorBarLimit, maxColorBarLimit)

    position=fig.add_axes([0.85, 0.15, 0.03, 0.70])
    cbar=plt.colorbar(cax=position)

    cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

    for t in cbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')   
        t.set_x(4.0)
    ###########################################################################
    #Add chain borders to the plot
    for i in range(len(selection_reorder)-1):
        beginningPoint = selection_reorder[i]/selection_reorder[-1]
        endingPoint = selection_reorder[i+1]/selection_reorder[-1]
        middlePoint = (float(beginningPoint)+float(endingPoint))/2.0
        if(i%2==0):
            #x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black'))

            #y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black'))  

        elif(i%2==1):
            #x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray'))
            
            #y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray')) 

        ax.annotate(myList[i], xy=(0, 1.04), xycoords='axes fraction', xytext=(middlePoint-0.015, 1.04), size=14, color='black')
        ax.annotate(myList[i], xy=(1.05, 0), xycoords='axes fraction', xytext=(1.05, middlePoint-0.015), rotation=90, size=14, color='black')
        #print(middlePoint)
    ###########################################################################
    #plt.tight_layout()       
    plt.savefig(out_file+'-overall-difference.png', bbox_inches='tight', dpi=200)
    #plt.show()

def projectCorrelationsOntoProteinVMD(ccMatrix, vmd_out_file, \
                                        selectedAtoms, valueFilter,\
                                        absoluteValues: bool, \
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
    #Calculate distance matrix
    dist_matrix=buildDistMatrix(selectedAtoms)
    
    #Plot the figure
    #print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    x=dist_matrix.flatten()
    y=ccMatrix.flatten()

    #print(len(y))

    distanceFilter = 0.5
    #valueFilter = 0.5
    #Write output in VMD format
    #Writing the output is very important for further analyses such as 
    #inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    #
    draw_string = "draw cylinder"+\
                " [lindex [[atomselect top \"chain {0:s} and resid {1:d} and name CA\"] get {{x y z}}] 0]"+\
                " [lindex [[atomselect top \"chain {2:s} and resid {3:d} and name CA\"] get {{x y z}}] 0]"+\
                " radius {4:.3f}\n"
    vdw_representation_string = "mol representation VDW 0.750000 50.000000\n"+\
                                "mol material Glossy\n"+\
                                "#mol color ColorID 7\n"+\
                                "mol selection \"chain {0:s} and resid {1:d} and name CA\"\n"+\
                                "mol addrep 0\n"
    DATA_FILE = open(vmd_out_file+'-general.tcl', 'w')


    DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
    DATA_FILE.write("mol modcolor 0 0 Chain\n")
    DATA_FILE.write("mol modmaterial 0 0 MetallicPastel\n")
    for i in range(0, len(ccMatrix)):
        for j in range(i+1, len(ccMatrix)):
            if(np.absolute(ccMatrix[i][j])>valueFilter):
                DATA_FILE.write(vdw_representation_string.\
                format(selectedAtoms.getChids()[i],\
                        selectedAtoms.getResnums()[i]))
                DATA_FILE.write(vdw_representation_string.\
                format(selectedAtoms.getChids()[j],\
                        selectedAtoms.getResnums()[j]))
                DATA_FILE.write(draw_string.format(selectedAtoms.getChids()[i],\
                        selectedAtoms.getResnums()[i],\
                        selectedAtoms.getChids()[j],\
                        selectedAtoms.getResnums()[j],\
                        #The radius of the connecting cylinder is proportional to the correlation value.
                        #However, it is necessary to multiply the radius with 0.5 to make it look better.
                        ccMatrix[i][j]*0.5))
    DATA_FILE.close()

    chains = Counter(selectedAtoms.getChids()).keys()

    plotChains = True
    if((len(chains)>1) & (plotChains)):
        #Inter-chain
        for chainI in chains:
            for chainJ in chains:
                if(chainI != chainJ):                
                    DATA_FILE = open(vmd_out_file+'-interchain-chains'+chainI+'-'+chainJ+'.tcl', 'w')
                    DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
                    DATA_FILE.write("mol modcolor 0 0 Chain\n")
                    DATA_FILE.write("mol modmaterial 0 0 MetallicPastel\n")
                    for i in range(0, len(ccMatrix)):
                        for j in range(i+1, len(ccMatrix)):
                            if(np.absolute(ccMatrix[i][j])>valueFilter):
                                if((selectedAtoms.getChids()[i] == chainI) and \
                                    (selectedAtoms.getChids()[j] == chainJ)):
                                    DATA_FILE.write(vdw_representation_string.\
                                    format(selectedAtoms.getChids()[i],\
                                            selectedAtoms.getResnums()[i]))
                                    DATA_FILE.write(vdw_representation_string.\
                                    format(selectedAtoms.getChids()[j],\
                                            selectedAtoms.getResnums()[j]))
                                    
                                    DATA_FILE.write(draw_string.\
                                    format(selectedAtoms.getChids()[i],\
                                            selectedAtoms.getResnums()[i],\
                                            selectedAtoms.getChids()[j],\
                                            selectedAtoms.getResnums()[j],\
                                            #The radius of the connecting cylinder is proportional to the correlation value.
                                            #However, it is necessary to multiply the radius with 0.5 to make it look better.
                                            ccMatrix[i][j]*0.5))
                    DATA_FILE.close()

        #Intra-chain
        for chain in chains:
            DATA_FILE = open(vmd_out_file+'-intrachain-chain'+chain+'.tcl', 'w')
            DATA_FILE.write("mol modstyle 0 0 NewCartoon 0.300000 50.000000 3.250000 0\n")
            DATA_FILE.write("mol modcolor 0 0 Chain\n")
            DATA_FILE.write("mol modmaterial 0 0 MetallicPastel\n")
            for i in range(0, len(ccMatrix)):
                for j in range(i+1, len(ccMatrix)):
                    if(np.absolute(ccMatrix[i][j])>valueFilter):
                        if((selectedAtoms.getChids()[i] == chain) and \
                           (selectedAtoms.getChids()[j] == chain)):
                            DATA_FILE.write(vdw_representation_string.\
                            format(selectedAtoms.getChids()[i],\
                                    selectedAtoms.getResnums()[i]))
                            DATA_FILE.write(vdw_representation_string.\
                            format(selectedAtoms.getChids()[j],\
                                            selectedAtoms.getResnums()[j]))    
                            DATA_FILE.write(draw_string.\
                            format(selectedAtoms.getChids()[i],\
                                    selectedAtoms.getResnums()[i],\
                                    selectedAtoms.getChids()[j],\
                                    selectedAtoms.getResnums()[j],\
                #The radius of the connecting cylinder is proportional to the correlation value.
                #However, it is necessary to multiply the radius with 0.5 to make it look better.
                                    ccMatrix[i][j]*0.5))
            DATA_FILE.close()
def correlationMapApp():
    print("@> Running 'Correlation Map App'")
    (inp_file, out_file, sel_type, pdb_file) = handle_arguments_correlationMaps()
    print("\n@> Input file   :", inp_file)
    print("@> PDB file     :", pdb_file)
    print("@> Data type    :", sel_type)    
    print("@> Output       :", out_file)

    ##########################################################################
    #Read PDB file 
    #TODO: This is the only place where I use Prody.
    #Maybe, I can replace it with a library that only parses 
    #PDB files. Prody does a lot more!
    selectedAtoms = parsePDB(pdb_file, subset='ca')

    ##########################################################################
    #Read data file and assign to a numpy array
    if(sel_type=="dcc"):
        ccMatrix=np.loadtxt(inp_file, dtype=float)
    elif(sel_type=="absdcc"):
        ccMatrix=np.absolute(np.loadtxt(inp_file, dtype=float))
    elif(sel_type=="lmi"):
        ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=True)
    else:
        print("Unknown matrix format!\n")
        sys.exit(-1)
    
    #Check the data type in the matrix.
    minCorrelationValue = np.min(ccMatrix)

    maxCorrelationValue = np.max(ccMatrix)

    if (minCorrelationValue<0.0):
        #Assume that it is an nDCC file
        minColorBarLimit = -1.0
    else:
        #Assume that it is an LMI file
        minColorBarLimit = 0.0

    if(maxCorrelationValue>1.0):
        print("This correlation map is not normalized!")
        #TODO: At this point, one can ask the user if s/he wants to normalize it!
        sys.exit(-1)
    ##########################################################################
    #Call overall correlation calculation
    maxColorBarLimit = 1.0
    overallCorrelationMap(ccMatrix, minColorBarLimit, maxColorBarLimit,\
                                    out_file, " ", selectedAtoms)

    plotDistributions = True
    if(plotDistributions):
        if(sel_type=="dcc"):
            distanceDistribution(ccMatrix, out_file, "nDCC", selectedAtoms, \
                                absoluteValues=False , writeAllOutput=True)

        elif(sel_type=="absdcc"):
            distanceDistribution(ccMatrix, out_file, "Abs(nDCC)", \
                selectedAtoms, absoluteValues=True , writeAllOutput=False)

        elif(sel_type=="lmi"):
            distanceDistribution(ccMatrix, out_file, "LMI", selectedAtoms, \
                                absoluteValues=True , writeAllOutput=True)

        else:
            print("Warning: Unknows correlation data.\n")
            print("         Correlations can be ndcc, absndcc, lmi!\n")

    ##########################################################################
    #Check number of chains. If there are multiple chains, plot inter and 
    #intra chain correlations
    chains = Counter(selectedAtoms.getChids()).keys()
    saveMatrix = False
    plotChains = True
    if((len(chains)>1) & (plotChains == True)):
        intraChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,\
                                        out_file, " ", selectedAtoms, saveMatrix)
        interChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,\
                                        out_file, " ", selectedAtoms, saveMatrix)

    #Here, we can filter some correlation values closer than a distance.
    #Typically, it is supposed to filter out the correlation within the 
    #same secondary structure etc. 
    filterByDistance = False
    if (filterByDistance):
        distanceValue=5.0
        ccMatrix=filterCorrelationMapByDistance(ccMatrix, out_file, " ", \
                                            selectedAtoms, distanceValue,\
                                            absoluteValues=False, \
                                            writeAllOutput=False)
    
    
    #Overall projection
    projectCorrelationsOntoProteinVMD(ccMatrix, out_file, 
                                        selectedAtoms, valueFilter=0.75,\
                                        absoluteValues=True, writeAllOutput=True)

    print("\n@> Program finished successfully!\n")
def usage_diffMaps():
    """
    Show how to use this program!
    """
    print("\nExample minimal usage:\n")
    print("python diffMap.py -i 4z90-cross-correlations.txt -j 4z91-cross-correlations.txt -p 4z90.pdb \n")
    print("Arguments: -i: The first file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)")
    print("           -j: The second file containing normalized dynamical cross correlations or LMI in matrix format. (Mandatory)")
    print("           -p: PDB file of the protein. (Mandatory)")
    print("           -q: A second PDB file for the other conformation if residues numbers are not same in two conformations. (Optional)")
    print("           -t: It can be ndcc, lmi or absndcc (absolute values of ndcc). Default value is ndcc (Optional)")
    print("           -o: This will be your output file. Output figures are in png format. (Optional)\n\n")

def handle_arguments_diffMaps():
    opts, args = getopt.getopt(sys.argv[1:],"hi:j:o:t:p:q:",["inp1=", "inp2=", "out=", "type=", "pdb=", "pdb2="])
    inp_file1 = None
    inp_file2 = None
    pdb_file1 = None
    pdb_file2 = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:j:o:t:p:q:",["inp1=", "inp2=", "out=", "type=", "pdb=", "pdb2="])
    except getopt.GetoptError:
        usage_diffMaps()
    for opt, arg in opts:
        if opt == '-h':
            usage_diffMaps()
            sys.exit(-1)
        elif opt in ("-i", "--inp1"):
            inp_file1 = arg
        elif opt in ("-j", "--inp2"):
            inp_file2 = arg    
        elif opt in ("-o", "--out"):
            out_file = arg
        elif opt in ("-t", "--type"):
            sel_type = arg
        elif opt in ("-p", "--pdb"):
            pdb_file1 = arg
        elif opt in ("-q", "--pdb2"):
            pdb_file2 = arg
        else:
            assert False, usage_diffMaps()

    #Input data matrix and PDB file are mandatory!
    if inp_file1==None or inp_file2==None or pdb_file1==None:
        usage_diffMaps()
        sys.exit(-1)

    #Assign a default name if the user forgets the output file prefix.
    if (out_file == None):
        out_file = "diff-map"

    #The user may prefer not to submit a title for the output.
    if (sel_type == None):
        sel_type = "ndcc"

    return (inp_file1, inp_file2, out_file, sel_type, pdb_file1, pdb_file2)

def diffMapApp():
    """
    This app helps to plot difference maps.
    It doesnt' assume that the the matrix sizes are equal in both maps.
    Moreover, if you provide a second pdb file, it can match residue 
    numbers and names of Calpha atoms in the pdb files. 
    Please, beware that the program does not do any sequence alignment.
    Therefore, if one of the proteins is contains a mutation or just a 
    relative of the first one, the app won't work. 
    """
    print("@> Running 'Difference Map App'")
    (inp_file1, inp_file2, out_file, sel_type, pdb_file1, pdb_file2) = handle_arguments_diffMaps()

    selectedAtomSet1 = parsePDB(pdb_file1, subset='ca')


    if (sel_type == 'lmi'):
        ccMatrix1 = convertLMIdata2Matrix(inp_file1, writeAllOutput=True)
        ccMatrix2 = convertLMIdata2Matrix(inp_file2, writeAllOutput=True)

        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif ((sel_type == 'absndcc')):
        ccMatrix1 = np.absolute(np.loadtxt(inp_file1, dtype=float))
        ccMatrix2 = np.absolute(np.loadtxt(inp_file2, dtype=float))
        minColorBarLimit = -1
        maxColorBarLimit = 1

    elif ((sel_type == 'ndcc')):
        ccMatrix1 = np.loadtxt(inp_file1, dtype=float)
        ccMatrix2 = np.loadtxt(inp_file2, dtype=float)
        minColorBarLimit = -2
        maxColorBarLimit = 2

    else:
        print("Error: Unknown matrix type!")
        print("What is the type of your correlation map: lmi, ndcc or absndcc?")
        sys.exit(-1)

    #One has to check if the lengths of two matrices match each other. 
    #Otherwise, there is a huge problem here. We have to match corresponding
    #residues and then do the subtraction. 
    if((len(ccMatrix1) == len(ccMatrix2)) and \
        (len(ccMatrix1[0]) == len(ccMatrix2[0]))):
        ccDiffMatrix = np.subtract(ccMatrix1, ccMatrix2)
        overallUniformDifferenceMap(ccMatrix1, ccMatrix2, minColorBarLimit, maxColorBarLimit, \
                                                                out_file, " ", selectedAtomSet1)
        #overallCorrelationMap(ccDiffMatrix, minColorBarLimit, maxColorBarLimit, \
        #                                                out_file, " ", selectedAtomSet1)
        ##########################################################################
        #Check number of chains. If there are multiple chains, plot inter and 
        #intra chain correlations
        chains = Counter(selectedAtomSet1.getChids()).keys()
        saveMatrix = False
        plotChains = True
        if((len(chains)>1) & (plotChains == True)):
            intraChainCorrelationMaps(ccDiffMatrix, minColorBarLimit, maxColorBarLimit,\
                                            out_file, " ", selectedAtomSet1, saveMatrix)
            interChainCorrelationMaps(ccDiffMatrix, minColorBarLimit, maxColorBarLimit,\
                                            out_file, " ", selectedAtomSet1, saveMatrix )
    else:
        print("@> Warning: Sizes of two matrices are not equal!")
        
        if(pdb_file2 == None):
            print("@> Warning: You have to specify at least two pdb files when")
            print("@>          matrix sizes are not identical!")
        else:
            selectedAtomSet2 = parsePDB(pdb_file2, subset='ca')
            overallNonUniformDifferenceMap(ccMatrix1, ccMatrix2, \
                                minColorBarLimit, maxColorBarLimit, \
                                out_file, " ", \
                                selectedAtomSet1, selectedAtomSet2)
def filterCorrelationMapByDistance(ccMatrix, out_file, title, \
                                    selectedAtoms, distanceValue,\
                                    absoluteValues: bool, writeAllOutput: bool):
    """
    If residues are closer to each other than a certain distance 
    (distanceValue), make these correlations zero. This filter can be useful 
    to get rid of short distance correlation for visualization purposes. 
    This function returns a filtered ccMatrix.
    """
    print("@> Filtering correlations lower than "+str(distanceValue)+" Angstrom")
    #Calculate distance matrix
    dist_matrix=buildDistMatrix(selectedAtoms)
    
    #print("@> Min. distance: {0:.2f} Angstrom.".format(np.min(dist_matrix)))
    #print("@> Max. distance: {0:.2f} Angstrom.".format(np.max(dist_matrix)))

    for i in range(0, len(ccMatrix)):
        for j in range(i+1, len(ccMatrix)):
            if(dist_matrix[i][j]<distanceValue):
                ccMatrix[i][j] = 0.0
                ccMatrix[j][i] = 0.0

    if(absoluteValues):
        dst_file = out_file+'-absolute-correlation-filtered'

    else:
        dst_file = out_file+'-correlation-filtered'
    
    #Write output
    #Writing the output is very important for further analyses such as 
    #inter-chain (inter-domain) or intra-chain (intra-domain) distributions etc.
    if(writeAllOutput):
        DATA_FILE = open(dst_file+'filtered.dat', 'w')
        for i in range(0, len(ccMatrix)):
            for j in range(i+1, len(ccMatrix)):
                DATA_FILE.write("{0:d}\t{1:s}\t{2:d}\t{3:s}\t{4:.3f}\t{5:.3f}\n".\
                                format(selectedAtoms.getResnums()[i],\
                                        selectedAtoms.getChids()[i],\
                                        selectedAtoms.getResnums()[j],\
                                        selectedAtoms.getChids()[j],\
                                        dist_matrix[i][j],\
                                        ccMatrix[i][j]))
        DATA_FILE.close()
    return ccMatrix

def networkAnalysis(ccMatrix, valueFilter, out_file, centrality, selectedAtoms):
    """
    This function calculates various network (graph) parameters of a protein.

    This function calculates some network centrality measures such as 
        -degree
        -betweenness
        -closeness
        -current flow betweenness
        -eigenvector.
    
    Parameters
    ----------
    ccMatrix: Numpy matrix
        It is a numpy matrix of typically nDCC, LMI or Generalized Correlations.
    valueFilter: float
        The ccMatrix values than the valueFilter will be ignored.
    out_file: string
        Prefix of the output file. According to the centralty measure, it will be 
        extended. 
    centrality: string
        It can have 'degree', 'betweenness', 'closeness' or 
        'current_flow'. 
    selectedAtoms: object
        This is a prody.parsePDB object of typically CA atoms of a protein.
    
    Returns
    -------
    Nothing
    """
    #Create your  graph
    dynNetwork = nx.Graph()
    

    n = selectedAtoms.numAtoms()

    #Add all CA atoms as nodes
    for i in range(n):
        dynNetwork.add_node(i)

    #Add all pairwise interactions greater than the valueFilter as edges.
    for i in range(n):
        for j in range(n):
            if(fabs(ccMatrix[i][j])>valueFilter):
                dynNetwork.add_edge(i, j, weight=-log(fabs(ccMatrix[i][j])))
                #dynNetwork.add_edge(i, j, weight=fabs(correlationArray[i][j]))

    ##########################Calculate degrees of all nodes
    if ((centrality == 'degree') or (centrality == 'all')):
        degreeResult = dynNetwork.degree(weight='weight')
        print("Degree calculation finished!")

        #open a file for degree
        degreeFile = open(out_file+"_degree_value_filter"+"{:.2f}".format(valueFilter)+'.dat', "w") 
        for i in range(n): 
        #    print(str(i)+" "+(str(dynNetwork.degree(i, weight='weight'))))
            degreeFile.write("{0:d}\t{1:.3f}\n".format((i)+1, degreeResult[i]))
        degreeFile.close()


    ##########################Calculate betweenness
    elif ((centrality == 'betweenness') or (centrality == 'all')):
        betweennessResult = nx.betweenness_centrality(dynNetwork, k=None, \
            normalized=True, weight='weight', endpoints=False, seed=None)

        #open a file for betweenness
        betweennessFile = open(out_file+"_betweenness_value_filter"+"{:.2f}".\
            format(valueFilter)+'.dat', "w") 
        print("Betweenness calculation finished!")
        for i in range(n): 
        #    print(str(i)+" "+(str(dynNetwork.betweenness(i, weight='weight'))))
            betweennessFile.write("{0:d}\t{1:.6f}\n".format((i)+1, betweennessResult[i]))
        betweennessFile.close()

    ##########################Calculate closeness
    elif ((centrality == 'closeness') or (centrality == 'all')):
        closenessResult = nx.closeness_centrality(dynNetwork, u=None, distance='weight')
        #open a file for closeness
        closenessFile = open(out_file+"_closeness_value_filter"+"{:.2f}".format(valueFilter)+'.dat', "w") 
        print("Closeness calculation finished!")
        for i in range(n): 
        #    print(str(i)+" "+(str(dynNetwork.closeness(i, weight='weight'))))
            closenessFile.write("{0:d}\t{1:.6f}\n".format((i)+1, closenessResult[i]))
        closenessFile.close()

    ##########################Calculate current_flow_betweenness
    elif ((centrality == 'current_flow') or (centrality == 'all')):
        current_flow_betweennessResult = nx.current_flow_betweenness_centrality(dynNetwork, normalized=True, weight='weight')

        #open a file for current_flow betweenness
        current_flow_betweennessFile = open(out_file+"_current_flow_betweenness_value_filter"+"{:.2f}".format(valueFilter)+'.dat', "w") 
        print("Current flow betweenness calculation finished!")
        for i in range(n): 
        #    print(str(i)+" "+(str(dynNetwork.betweenness(i, weight='weight'))))
            current_flow_betweennessFile.write("{0:d}\t{1:.6f}\n".format((i)+1, current_flow_betweennessResult[i]))
        current_flow_betweennessFile.close()

if __name__ == "__main__":
        #TODO:
    # There are a bunch of things one can add to this script:
    # 1-Plot nDCC maps or normalized linear mutual information maps!: Done!
    # 2-Project (high) correlations onto PDB structure.
    #   a) as a Pymol script output
    #   b) as a VMD script output: Done!
    # 3-Project secondary structures on x and y axes of a correlation map.
    # 4-Difference maps: Done!
    # 5-Combining two correlation plots as upper triangle and lower triangle. 
    # 6-Filter correlations lower than a certain (absolute) value.: Done!
    # 7-Filter correlations for residues that are very close.: Done!
    # 8
    print("\n\n|------------------------------Correlation Plus------------------------------|")
    print("|                                                                            |")
    print("|   A set of utility programs to plot and analyze protein correlation maps.  |")
    print("|               Copyright (c) 2019-2020 Mustafa Tekpinar                     |")
    print("|                       Email: tekpinar@buffalo.edu                          |")
    print("|                          Licence: MIT License                              |")
    print("|--------------------------------------------------------------------------- |\n\n")
    correlationMapApp()

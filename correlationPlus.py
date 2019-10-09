#!/usr/bin/env python
"""
Program Name: Cross-correlation Plotting Program (I'll find a fancy name later!)
Author      : Mustafa TEKPINAR
Email       : tekpinar@buffalo.edu
Copyright   : Mustafa Tekpinar - 2019
License     : BSD

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
from collections import Counter

def usage():
    """
    Show how to use the small program!
    """
    print("\nExample usage:\n")
    print("python nDCCPotGeneral.py -i 4z90_cross-correlations.txt -s ' ' -p 4z90.pdb -o 4z90_cross-correlations\n")
    print("Arguments: -i: A file containing normalized dynamical cross correlations in matrix format.")
    print("           -s: A string for the title of the plot.")
    print("           -p: PDB file of the protein.")
    print("           -o: This will be your output file. Output is in png format.")

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
def handle_arguments():
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:s:p:",["inp=", "out=", "sel=", "pdb="])
    inp_file = None
    pdb_file = None
    out_file = None
    sel_type = None

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:s:p:",["inp=", "out=", "sel=", "pdb="])
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
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
            assert False, usage()

    if inp_file==None or out_file==None:
        usage()
        sys.exit(-1)
    return (inp_file, out_file, sel_type, pdb_file)

def overall_nDCC_map(ccMatrix, out_file, sel_type, selectedAtoms):
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
    plt.title(sel_type, y=1.08)

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
    #print(np.argmin(ccMatrix_sub, 0))

    #print("Column Index of min value")
    #print(np.argmin(ccMatrix_sub, 1))

    #Set colorbar features here!  
    jet=plt.get_cmap('jet') 
    djet = cmap_discretize(jet, 8)

    plt.imshow(np.matrix(ccMatrix), cmap=djet)
    plt.clim(-1.0, 1.0)

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
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray'))  

        elif(i%2==1):
            #x axis
            ax.annotate('', xy=(beginningPoint, 1.03), xycoords='axes fraction', \
                        xytext=(endingPoint, 1.03), arrowprops=dict(linewidth = 2., arrowstyle="-", color='gray'))
            
            #y axis
            ax.annotate('', xy=(1.04, beginningPoint), xycoords='axes fraction', \
                        xytext=(1.04, endingPoint), arrowprops=dict(linewidth = 2., arrowstyle="-", color='black')) 

        ax.annotate(myList[i], xy=(0, 1.04), xycoords='axes fraction', xytext=(middlePoint-0.015, 1.04), size=14, color='black')
        ax.annotate(myList[i], xy=(1.05, 0), xycoords='axes fraction', xytext=(1.05, middlePoint-0.015), rotation=90, size=14, color='black')
        #print(middlePoint)

    plt.tight_layout()       
    plt.savefig(out_file+'.png', dpi=200)
    #plt.show()

def intrachain_nDCC_maps(ccMatrix, out_file, sel_type, selectedAtoms):
    """
    Plot intra-chain correlations if there are at least two chains!
    """

    myList = list(Counter(selectedAtoms.getChids()).keys())
    print(myList)

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

    print(selection_reorder)

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
        plt.clim(-1.0, 1.0)

        position=fig.add_axes([0.85, 0.15, 0.03, 0.70])
        
        cbar=plt.colorbar(cax=position)
        cbar.set_ticks([-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])

        for t in cbar.ax.get_yticklabels():
            t.set_horizontalalignment('right')   
            t.set_x(4.0)
        
        plt.savefig(out_file+'-chain'+myList[j]+'.png', dpi=200)

        #plt.show()
        plt.close('all')

def interchain_nDCC_maps(ccMatrix, out_file, sel_type, selectedAtoms):
    """
    Plot inter-chain correlations if there are at least two chains!
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    myList = list(Counter(selectedAtoms.getChids()).keys())
    print(myList)
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

    print(selection_reorder)

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
            plt.axis([0, (selection_reorder[l+1]-selection_reorder[l]), 0, (selection_reorder[k+1]-selection_reorder[k])])

            sub_nDCC_matrix = ccMatrix[(selection_reorder[k]) : (selection_reorder[k+1]), (selection_reorder[l]) : (selection_reorder[l+1])]
            plt.imshow(np.matrix(sub_nDCC_matrix), cmap=djet)
            plt.clim(-1.0, 1.0)

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

if __name__ == "__main__":
    #TODO:
    # There are a bunch of things one can do with this script:
    # 1-Plot nDCC maps or normalized linear mutual information maps!
    # 2-Project (high) correlations onto PDB structure.
    #   a) as a pymol script output
    #   b) as a VMD script output
    # 3-Project secondary structures on x and y axes.

    (inp_file, out_file, sel_type, pdb_file) = handle_arguments()
    print("\n@> Input file   :", inp_file)
    print("@> Pdb file     :", pdb_file)
    print("@> Output       :", out_file)
    print("@> Your title   :", sel_type)

    ##########################################################################
    #Read PDB file 
    #TODO: This is the only place where I use Prody.
    #Maybe, I can replace it with a library that only parses 
    #PDB files. Prody does a lot more!
    selectedAtoms = parsePDB(pdb_file, subset='ca')

    ##########################################################################
    #Read data file and assign to a numpy array
    ccMatrix=np.loadtxt(inp_file, dtype=float)
    
    ##########################################################################
    #Call overall nDCC calculation
    overall_nDCC_map(ccMatrix, out_file, sel_type, selectedAtoms)

    ##########################################################################
    #Check number of chains. If there are multiple chains, plot inter and 
    #intra chain correlations
    chains = Counter(selectedAtoms.getChids()).keys()
    if(len(chains)>1):
        intrachain_nDCC_maps(ccMatrix, out_file, sel_type, selectedAtoms)
        interchain_nDCC_maps(ccMatrix, out_file, sel_type, selectedAtoms)

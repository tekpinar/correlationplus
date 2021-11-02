###############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   correlations maps of proteins.                 #
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


"""
    One of the most important advantages of correlationplus compared to webservers
    performing similar tasks is its flexibility provided by its Application 
    Programming Interface (API). 
    Here, we provide two examples of API usage. 
"""
#############EXAMPLE 1:
# The first example shows how you can modify the output visualizations according 
# to your own taste and needs. Here, we utilize a high-quality visualization
# library called seaborn. 

from correlationplus.calculate import calcMD_LMI, calcENM_LMI
import prody
import matplotlib.pylab as plt
import seaborn as sbrn

# Parse the necessary files that are already in the examples folder 
protein = "6lu7_dimer_with_N3_protein_sim1_ca.pdb"
selectedAtoms=prody.parsePDB(protein, subset='ca')

# Call the calculate function in correlationplus
LMImatrix = calcENM_LMI(selectedAtoms, cut_off=15, saveMatrix=False)

# Use seaborn for a continuous heatmap instead of the default 
# discrete heatmap implemented in the correlationplus.  
sbrn.heatmap(LMImatrix, vmin=0.0, vmax=1.0, cmap='viridis')
plt.show()

#############EXAMPLE 2:
# In the 2nd example, we will calculate normalized linear mutual informatio (LMI)
# as a moving average to see how conformational changes may alter the dynamical
# correlations. The produced images will be converted to a gif file to help the analysis. 
# We will use the trajectory file given in the examples folder of the 
# correlationplus repository.  
from correlationplus.visualize import overallCorrelationMap
trajectory = "6lu7_dimer_with_N3_protein_sim1_ca_short.trr"
windowSize = 100
filenames = []

# Compute the LMI maps with a windowSize=100
# We skip every 10 steps to make the calculations faster.
for i in range(0, 501-windowSize, 10):
    outfile = "LMI_frame_"+str(i)
    filenames.append(outfile+"-overall.png")
    LMImatrix = calcMD_LMI(protein, trajectory, \
                    startingFrame=i, endingFrame=(i+100), \
                    normalized=True, alignTrajectory=True, \
                    atomSelection='protein and name CA', \
                    saveMatrix=False) 
    overallCorrelationMap(LMImatrix, 0.0, 1.0, outfile, " ", selectedAtoms)

# Make a gif file from the output images
import imageio 
images = [] 
for filename in filenames: 
    images.append(imageio.imread(filename))
imageio.mimsave('movie.gif', images)

# Clean all unnecessary files
import os
for filename in filenames: 
    os.remove(filename)

# As shown the examples above, the API may allow more systematic analyses quite 
# efficiently. Moreover, one can integrate quickly various Python libraries both to test new ideas or 
# improve visualizations according to the needs. The Python API of the correlationplus brings with itself
# a significant versatily compared to the webservers that perform similar tasks. 

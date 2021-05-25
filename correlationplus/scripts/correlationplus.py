###############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                    dynamical correlations maps of proteins.                 #
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

# cannot use relative imports
# as it does not work in tests/run_test
# ImportError: attempted relative import with no known parent package
from correlationplus.scripts.calculate import calculateApp
from correlationplus.scripts.visualize import visualizemapApp
from correlationplus.scripts.diffMap import diffMapApp
from correlationplus.scripts.analyze import centralityAnalysisApp
from correlationplus.scripts.paths import pathAnalysisApp

#TODO:
# There are a bunch of things one can add to this script:
# 1-Plot nDCC maps or normalized linear mutual information maps!: Done!
# 2-Project (high) correlations onto PDB structure.
#   a) as a Pymol script output: Done
#   b) as a VMD script output: Done!
# 3-Project secondary structures on x and y axes of a correlation map.
# 4-Difference maps: Done!
# 5-Combining two correlation plots as upper triangle and lower triangle.
# 6-Filter correlations lower than a certain (absolute) value.: Done!
# 7-Filter correlations for residues that are very close.: Done!
# 8-Add centrality calculations: Done
# 9-Add centrality visualizations: Done
# 10-Time lagged correlation analysis for MD trajectories:
# 11-Time shifted correlation analysis. For example, shift 100 frames each time
#    to do centrality or community analysis.
# 12-Better pruning functions to clean correlation maps.
# 13-A new utilities module to add, subtract, average, triangulate the correlation 
#    maps. Some of the currently available functionalities can be moved there.
# 14-A better way to run tests more systematically.
# 15-Move to argparse instead of getopt.  


def usage_main():
    """
    Show how to use this program!
    """
    print("""
Example usage:

correlationplus -h

CorrelationPlus contains four apps:
 - calculate
 - visualize
 - analyze
 - diffMap

You can get more information about each individual app as follows:

correlationplus analyze -h
""")


def main():

    print("""

|------------------------------Correlation Plus------------------------------|
|                                                                            |
|        A Python package to calculate, visualize and analyze protein        |
|                           correlation maps.                                |
|               Copyright (C) Mustafa Tekpinar, 2017-2018                    |
|                   Copyright (C) CNRS-UMR3528, 2019                         |
|             Copyright (C) Institut Pasteur Paris, 2020-2021                |
|                         Author: Mustafa Tekpinar                           |
|                       Email: tekpinar@buffalo.edu                          |
|                           Licence: GNU LGPL V3                             |
|--------------------------------------------------------------------------- |

""")

    if len(sys.argv) > 1:
        if sys.argv[1] == "calculate":
            calculateApp()
        elif sys.argv[1] == "visualize":
            visualizemapApp()
        elif sys.argv[1] == "analyze":
            centralityAnalysisApp()
        elif sys.argv[1] == "paths":
            pathAnalysisApp()
        elif sys.argv[1] == "diffMap":
            diffMapApp()
        elif sys.argv[1] == "-h" or sys.argv[1] == "--help":
            usage_main()
        else:
            usage_main()
            sys.exit(-1)
    else:
        usage_main()
        sys.exit(-1)


if __name__ == "__main__":
    main()

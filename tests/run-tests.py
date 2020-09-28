###############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                    dynamical correlations maps of proteins.                 #
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

import os

"""
Unfortunately, this is not the most systematic way of writing unit tests.
But I needed a quick and dirty way of checking integrity of the apps.
So, I will do it here by checking diffs of output images.
"""


def runTests():
    # Test correlationMapApp for nDCC maps
    prefix = os.path.normpath(os.path.join(__file__, '..', '..'))

    # Test calculate module to produce nDCC maps with GNM.
    os.system(f"correlationplus calculate "
              f"-p {prefix}/examples/6fl9_centeredOrientedAligned2Z.pdb "
               "-m GNM -o gnm-ndcc.dat")

    # Test calculate module to produce nDCC maps with ANM.
    os.system(f"correlationplus calculate "
              f"-p {prefix}/examples/6fl9_centeredOrientedAligned2Z.pdb "
               "-m ANM -o anm-ndcc.dat")


    # Test visualizemapApp for absolute nDCC maps
    os.system(f"correlationplus visualizemap "
              f"-i {prefix}/examples/6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt "
              f"-p {prefix}/examples/6fl9_centeredOrientedAligned2Z.pdb "
               "-t absndcc")

    # Test visualizemapApp for LMI maps
    os.system(f"correlationplus visualizemap "
              f"-i {prefix}/examples/6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat "
              f"-p {prefix}/examples/6fl9_centeredOrientedAligned2Z.pdb "
              "-t lmi")

 
    # Test diffMapApp for LMI maps
    os.system(f"correlationplus diffMap "
              f"-i {prefix}/examples/6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat "
              f"-j {prefix}/examples/zacharias_rc15_scalCoeff15_100_modes_lmi.dat "
              f"-p {prefix}/examples/6fl9_centeredOrientedAligned2Z.pdb "
              "-t lmi")

    # Test centralityAnalysisApp for LMI maps
    os.system(f"correlationplus centralityAnalysis "
              f"-i {prefix}/examples/6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat "
              f"-p {prefix}/examples/6fl9_centeredOrientedAligned2Z.pdb "
              "-t lmi")

if __name__ == "__main__":
    import sys
    COR_PLUS_HOME = os.path.abspath(os.path.join(__file__, '..', '..'))
    old_path = sys.path
    if COR_PLUS_HOME not in sys.path:
        sys.path.insert(0, COR_PLUS_HOME)
    runTests()
    sys.path = old_path

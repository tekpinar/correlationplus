###############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                     correlation maps of proteins.                           #
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

import os

"""
Unfortunately, this is not the most systematic way of writing unit tests.
But I needed a quick and dirty way of checking integrity of the apps.
So, I will do it here by checking diffs of output images.
"""


def runTests():
    # Test correlationMapApp for nDCC maps
    prefix = os.path.normpath(os.path.join(__file__, '..', '..'))
    # Test calculate module to produce nDCC map with GNM.
    command = f"correlationplus calculate -p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb  " \
              f"-m GNM -o ndcc-6lu7-gnm.dat"
    print(f"\n@> Testing the nDCC calculation from GNM with the following command:\n\n{command}\n\n")
    os.system(command)

    # Test calculate module to produce nDCC map with ANM.
    command = f"correlationplus calculate -p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb  " \
              f"-m ANM -o ndcc-6lu7-anm.dat"
    print(f"\n@> Testing the nDCC calculation from ANM with the following command:\n\n{command}\n\n")
    os.system(command)

    # Test calculate module to produce nDCC (normalized dynamical cross-correlation) map from an MD trajectory.
    command = f"correlationplus calculate "\
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "\
              f"-f {prefix}/examples/6lu7_dimer_with_N3_protein_sim1_short.trr -o ndcc-6lu7-md.dat"
    print(f"\n@> Testing the nDCC calculation from an MD trajectory with the following command:\n\n{command}\n\n")
    os.system(command)

    # Test calculate module to produce DCC (dynamical cross-correlation) map from an MD trajectory.
    # Please note that DCC values are not normalized here!
    command = f"correlationplus calculate "\
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "\
              f"-f {prefix}/examples/6lu7_dimer_with_N3_protein_sim1_short.trr -o dcc-6lu7-md.dat -t dcc"
    print(f"\n@> Testing the DCC calculation from an MD trajectory with the following command:\n\n{command}\n\n")
    os.system(command)

    # Test calculate module to produce a normalized LMI map with ANM modes.
    os.system(f"correlationplus calculate "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
               "-t nlmi -o nlmi-6lu7-anm.dat")

    # Test calculate module to produce a NON-normalized LMI map with ANM modes.
    os.system(f"correlationplus calculate "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
               "-t lmi -o lmi-6lu7-anm.dat")

    # Test calculate module to produce NON-normalized LMI map from an MD trajectory.
    os.system(f"correlationplus calculate "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              f"-f {prefix}/examples/6lu7_dimer_with_N3_protein_sim1_short.trr -t lmi -o lmi-6lu7-md.dat")

    # Test calculate module to produce normalized LMI map from an MD trajectory.
    os.system(f"correlationplus calculate "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              f"-f {prefix}/examples/6lu7_dimer_with_N3_protein_sim1_short.trr -t nlmi -o nlmi-6lu7-md.dat")

    # Test visualizemapApp for absolute nDCC maps
    os.system(f"correlationplus visualize "
              f"-i {prefix}/examples/ndcc-6lu7-anm.dat "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb  "
              "-t absndcc -v 0.75")

    # Test visualizemapApp for LMI maps
    os.system(f"correlationplus visualize "
              f"-i {prefix}/examples/6lu7_dimer_with_N3_protein_sim1-lmi.dat "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              "-t nlmi -v 0.75 ")
 
    # Test diffMapApp for nLMI maps
    os.system(f"correlationplus diffMap "
              f"-i {prefix}/examples/6lu7_dimer_with_N3_protein_sim1-lmi.dat  "
              f"-j {prefix}/examples/6lu7_dimer_no_ligand_protein_sim1-lmi.dat  "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              "-t nlmi")

    # Test centralityAnalysisApp for LMI maps
    os.system(f"correlationplus analyze "
              f"-i {prefix}/examples/6lu7_dimer_with_N3_protein_sim1-lmi.dat "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              "-t nlmi")

    # Test pathAnalysisApp for ndcc maps
    os.system(f"correlationplus paths "
              f"-i {prefix}/examples/ndcc-6lu7-anm.dat "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              "-b A41 -e B41")

    # Test pathAnalysisApp for LMI maps
    os.system(f"correlationplus paths "
              f"-i {prefix}/examples/6lu7_dimer_with_N3_protein_sim1-lmi.dat "
              f"-p {prefix}/examples/6lu7_dimer_with_N3_protein_sim1.pdb "
              "-t nlmi -b A41 -e B41")

    # Test to visualize elasticity graphs (FitNMA) obtained from
    # https://onlinelibrary.wiley.com/doi/10.1002/jcc.26701?af=R
    os.system(f"correlationplus visualize -i {prefix}/examples/crn.enm -p {prefix}/examples/1crn.pdb -t eg -r 0.05")

    # Test to analyze elasticity graphs (FitNMA) obtained from
    # https://onlinelibrary.wiley.com/doi/10.1002/jcc.26701?af=R
    os.system(f"correlationplus analyze -i {prefix}/examples/crn.enm -p {prefix}/examples/1crn.pdb -t eg -v 0.1 -d 30")

    # Test to visualize paths in elasticity graphs (FitNMA) obtained from
    # https://onlinelibrary.wiley.com/doi/10.1002/jcc.26701?af=R
    os.system(f"correlationplus paths -i {prefix}/examples/crn.enm -p {prefix}/examples/1crn.pdb -t eg -v 0.1 -d 30 -b A10 -e A41")
    
    
if __name__ == "__main__":
    import sys
    COR_PLUS_HOME = os.path.abspath(os.path.join(__file__, '..', '..'))
    old_path = sys.path
    if COR_PLUS_HOME not in sys.path:
        sys.path.insert(0, COR_PLUS_HOME)
    runTests()
    sys.path = old_path

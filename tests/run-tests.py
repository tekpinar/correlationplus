import os
import sys
"""
Unfortunately, this is not the most systematic way of writing unit tests.
But I needed a quick and dirty way of checking integrity of the apps.
So, I will do it here by checking diffs of output images.
"""

def runTests():
    #Test correlationMapApp for nDCC maps
    os.system("python ../src/correlationPlus.py mapAnalysisApp"\
        +" -i ../examples/6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt"\
        +" -p ../examples/6fl9_centeredOrientedAligned2Z.pdb"\
        +" -s absdcc")

    #Test correlationMapApp for absolute nDCC maps
    #os.system()

    #Test correlationMapApp for LMI maps
    #os.system()

    #Test diffMapApp for LMI maps
    os.system("python ../src/diffMap.py diffMapApp"+\
        " -i ../examples/6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat"+\
        " -j ../examples/zacharias_rc15_scalCoeff15_100_modes_lmi.dat"+\
        " -p ../examples/6fl9_centeredOrientedAligned2Z.pdb"+\
        " -t lmi")

if __name__ == "__main__":
    runTests()

    

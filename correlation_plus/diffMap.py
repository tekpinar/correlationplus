"""
This module helps to plot difference maps.
It doesnt' assume that the the matrix sizes are equal in both maps.
Moreover, if you provide a second pdb file, it can match residue 
numbers and names of Calpha atoms in the pdb files. 
Please, beware that the program does not do any sequence alignment.
Therefore, if one of the proteins is contains a mutation or just a 
relative of the first one, the app won't work. 
"""

# from collections import Counter, OrderedDict
import correlationPlus as cP
# import numpy as np
# from prody import parsePDB
# import matplotlib.pyplot as plt

if __name__ == "__main__":
    cP.diffMapApp()
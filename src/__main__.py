#!/usr/bin/env python
"""
Program Name: Cross-correlation Plotting Program (I'll find a fancy name later!)
Author      : Mustafa TEKPINAR
Email       : tekpinar@buffalo.edu
Copyright   : Mustafa Tekpinar - 2019
License     : MIT License

Purpose     : This is a small program to automatize plotting of normalized 
dynamical cross-correlations obtained from molecular dynamics simulations or 
elastic network models. This script can be useful if you have multiple 
chains in a structure and you want to see intra-chain and inter-chain 
correlations more clearly. I just didn't like the way current programs are 
doing it and I wrote something for myself. I hope it may help the others also!
"""

from correlationPlus import *

if __name__ == "__main__":
    main()

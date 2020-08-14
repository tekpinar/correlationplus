###############################################################################
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

"""
Program Name: Cross-correlation Plotting Program (I'll find a fancy name later!)
Author      : Mustafa TEKPINAR
Email       : tekpinar@buffalo.edu

Purpose     : This is a small program to automatize plotting of normalized
dynamical cross-correlations obtained from molecular dynamics simulations or
elastic network models. This script can be useful if you have multiple
chains in a structure and you want to see intra-chain and inter-chain
correlations more clearly. I just didn't like the way current programs are
doing it and I wrote something for myself. I hope it may help the others also!
"""

__all__ = ['mapAnalysis', 'diffMap', 'centralityAnalysis']

__version__ = '0.0.5'


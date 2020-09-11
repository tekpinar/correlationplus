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

from setuptools import setup, find_packages

from correlationPlus import __version__ as cp_vers

setup(name='correlationPlus',
      version=cp_vers,
      description=" Python module to plot dynamical correlations maps for proteins.",
      long_description=open('README.md').read(),
      long_description_content_type="text/markdown",
      author="Mustafa Tekpinar",
      author_email="tekpinar@buffalo.edu",
      url="https://github.com/tekpinar/correlationPlus",
      download_url="https://github.com/tekpinar/correlationPlus",
      license="LGPL",
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
          ],
      python_requires='>=3.6',
      install_requires=open("requirements.txt").read().split(),
      # zip_safe=False,
      packages=[p for p in find_packages() if p != 'tests'],
      # file where some variables must be fixed by install
      entry_points={
          'console_scripts': [
              'correlationPlus=correlationPlus.scripts.correlationPlus:main'
          ]
      }
      )

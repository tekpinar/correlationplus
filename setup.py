from setuptools import setup, find_packages

from correlation_plus import __version__ as cp_vers

setup(name='correlationPlus',
      version=cp_vers,
      description=" Python module to plot dynamical correlations maps for proteins.",
      long_description=open('README').read(),
      author="Mustafa Tekpinar",
      author_email="tekpinar@buffalo.edu",
      url="",
      download_url='',
      license="MIT License",
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
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
              'mapAnalysis=correlation_plus.scripts.mapAnalysis:mapAnalysisApp',
              'diffMap=correlation_plus.scripts.diffMap:diffMapApp'
          ]
      }
      )

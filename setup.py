#from distutils.core import setup
from setuptools import setup

#This line was after license. 
##packages=['correlationPlus'],
setup(name='correlationPlus', 
      package_dir = {'correlationPlus': 'src'},
      entry_points={
        'console_scripts': [
            'correlationPlus = correlationPlus.__main__:main',
        ],
      },
#      scripts=['src/correlationPlus.py'],
      version='0.0.5', 
      description='A Python module to plot dynamical correlations for proteins.',
      author='Mustafa Tekpinar',
      author_email='tekpinar@buffalo.edu',
      license='MIT License',
      url='https://github.com/tekpinar/correlationPlus',
      install_requires=['numpy', 'prody', 'matplotlib'],
      )


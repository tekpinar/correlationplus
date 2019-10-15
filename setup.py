from distutils.core import setup

setup(name='correlationPlus', 
      package_dir = {'correlationPlus': 'src'},
      version='0.0.4', 
      description='A Python module to plot dynamical correlations for proteins.',
      author='Mustafa Tekpinar',
      author_email='tekpinar@buffalo.edu',
      license='MIT License',
      url='https://github.com/tekpinar/correlationPlus',
      packages=['correlationPlus'],
      install_requires=['numpy', 'prody', 'matplotlib', 'getopt'],
      )


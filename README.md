# correlationplus

A Python package to calculate, visualize and analyze correlation maps of proteins.

correlationplus contains four main scripts that you can use to calculate, visualize
and analyze correlation maps of proteins. 
These correlations can be dynamical cross-correlations, linear mutual
information, sequence based covariation/coevolution or any other pairwise coupling metric. 
You can use elastic network models or your molecular dynamics trajectories for calculation 
of dynamical correlations.  

## Video Tutorials
* Installation: https://www.youtube.com/watch?v=Fc_xpnrbbWU
* Calculate Module: https://www.youtube.com/watch?v=HG8Mzmy1IOM
* Visualize Module: https://www.youtube.com/watch?v=HgUVAV1unXs
* Analyze Module: https://www.youtube.com/watch?v=04BwJDauOn4
* Sequence Conservation Analysis with CoeViz and correlationplus (Part 1): https://www.youtube.com/watch?v=ieu91glEd0s
* Sequence Conservation Analysis with CoeViz and correlationplus (Part 2): https://www.youtube.com/watch?v=BuxQNFoid1A

## Installation

We recommend one of the installation methods for regular users:

### with conda
```bash
conda install -c bioconda correlationplus

```

### with pip

The pip version required by some dependencies is >= 21.0.1, which is not the pip version bundle with python 3.(6,7,8)
So, you have to update pip before installing *correlationplus*. Otherwise, you will have trouble during *MDAnalysis* dependency installation.
For this reason, we **strongly** encourage you to install correlationplus in a [virtualenv](https://virtualenv.pypa.io/en/latest/).

```bash
python3 -m venv correlationplus
cd correlationplus
source bin/activate
python3 -m pip install -U pip
python3 -m pip install correlationplus
```

if you want to install it without using a virtualenv
and encounter an error related to llvmlite, you can
solve it as follows:
```bash
python3 -m pip install llvmlite --ignore-installed
python3 -m pip install correlationplus
```

or if you do not have administration rights
```bash
python3 -m pip install --user llvmlite --ignore-installed
python3 -m pip install --user correlationplus
```

Most of the time, at least one these methods will be sufficient for the installation.
However, if these two methods didn't work for any reason, you can take a look 
to 'Advanced Installation' instructions at
https://correlationplus.readthedocs.io/en/latest/installation.html#advanced-installation.


## Quickstart
There are four main scrips: 
* calculate
* visualize
* analyze
* paths

You can find details of each script with the following commands:

```bash
correlationplus calculate -h
correlationplus visualize -h
correlationplus analyze -h
correlationplus paths -h
```

Go to our [readthedocs](https://correlationplus.readthedocs.io/en/latest/quickstart.html) for 
detailed usage examples for each script.

## Licensing

*correlationplus* is developed and released under [GNU Lesser GPL Licence](https://www.gnu.org/licenses/lgpl-3.0.en.html). 
Please read to the **COPYING** and **COPYING.LESSER** files to know more. 

# correlationplus

A Python package to calculate, visualize and analyze dynamical correlations of proteins.

correlationplus contains three main scripts that you can use to calculate, visualize
and analyze dynamical correlations of proteins. 
These correlations can be dynamical cross-correlations, linear mutual
information. You can use elastic network models or your
your molecular dynamics trajectories for calculation, visualizations and analyses. 

## Installation

We recommend one of the installation methods for regular users:

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

if you still want to install it without using a virtualenv
or you hate to upgrade pip version to a version >=21.0.1,
install numpy==1.16
and the install correlationplus

We recommend to use pip as follows:
```bash
python3 -m pip install numpy==1.16.1
python3 -m pip install correlationplus
```

or if you do not have administration rights
```bash
python3 -m pip install --user numpy==1.16.1
python3 -m pip install --user correlationplus
```


### with conda
```bash
conda install -c bioconda correlationplus

```

Most of the time, at least one these methods will be sufficient for the installation.
However, if these two methods didn't work for any reason, you can take a look 
to 'Advanced Installation' instructions at
https://correlationplus.readthedocs.io/en/latest/installation.html#advanced-installation.


## Quickstart
There are three main scrips: 
* calculate
* visualize
* analyze

You can find details of each script with the following commands:

```bash
correlationplus calculate -h
correlationplus visualize -h
correlationplus analyze -h
```

Go to our [readthedocs](https://correlationplus.readthedocs.io/en/latest/quickstart.html) for 
detailed usage examples for each script.

## Licensing

*correlationplus* is developed and released under [GNU Lesser GPL Licence](https://www.gnu.org/licenses/lgpl-3.0.en.html). 
Please read to the **COPYING** and **COPYING.LESSER** files to know more. 

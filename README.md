# correlationplus

A Python package to calculate, visualize and analyze dynamical correlations of proteins.

correlationplus contains three main scripts that you can use to calculate, visualize
and analyze dynamical correlations of proteins. 
These correlations can be dynamical cross-correlations, linear mutual
information or generalized correlations. You can use elastic network models or your
your molecular dynamics trajectories for these visualizations and analyses. 

## Installation

We recommend on of the installation methods for regular users:

### with pip

We recommend to use pip as follows:
```bash
pip install correlationplus
```

or if you do not have administration rights
```bash
pip install --user correlationplus
```

If you prefer to use a virtualenv
```bash
python3 -m venv correlationplus
cd correlationplus
source bin/activate
pip install correlationplus
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

You can find details of each script below.

```bash
correlationplus calculate -h
correlationplus visualize -h
correlationplus analyze -h
```

Go to our [readthedocs](https://correlationplus.readthedocs.io/en/latest/quickstart.html) for detailed usage examples.

## Licensing

*correlationplus* is developed and released under [GNU Lesser GPL Licence](https://www.gnu.org/licenses/lgpl-3.0.en.html). 
Please read to the **COPYING** and **COPYING.LESSER** files to know more. 
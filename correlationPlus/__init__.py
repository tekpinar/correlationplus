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


import pkgutil
import inspect
from pathlib import Path
from importlib import import_module


def _import_public_function():
    """
    import all function in al modules (except itself and scripts) from this package
    and inject them in the package namespace
    so this function can be used like this:
        import correlation_plus
        correlation_plus.findCommonCorePDB(...

    or ::

        from correlation_plus import *
        findCommonCorePDB(...
    """
    for (_, name, _) in pkgutil.iter_modules([Path(__file__).parent]):
        if not (name.startswith('__') or name == 'scripts'):
            imported_module = import_module('.' + name, package=__name__)
            pub_funcs = [(f_name, f) for f_name, f in inspect.getmembers(imported_module, inspect.isfunction)
                         if not f_name.startswith('_')]
            for f_name, func in pub_funcs:
                globals()[f_name] = func


__version__ = '0.0.5'

_import_public_function()

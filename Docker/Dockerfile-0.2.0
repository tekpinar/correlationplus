##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   dynamical correlations maps of proteins.                  #
# Authors: Mustafa Tekpinar                                                   #
# Copyright Mustafa Tekpinar 2017-2018                                        #
# Copyright CNRS-UMR3528, 2019                                                #
# Copyright Institut Pasteur Paris, 2020-2021                                 #
#                                                                             #
# This file is part of correlationplus.                                       #
#                                                                             #
# correlationplus is free software: you can redistribute it and/or modify     #
# it under the terms of the GNU Lesser General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# correlationplus is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU LESSER General Public License for more details.                         #
#                                                                             #
# You should have received a copy of the GNU Lesser General Public License    #
# along with correlationplus.  If not, see <https://www.gnu.org/licenses/>.   #
###############################################################################
FROM python:3.8.5


LABEL org.label-schema.vcs-url="https://github.com/tekpinar/correlationplus"
LABEL org.label-schema.version="0.2.0"
LABEL org.label-schema.description="A Python package to calculate, visualize and analyze dynamical correlations of proteins."
LABEL org.label-schema.docker.cmd="docker run -v ~:/home/correlationplus correlationplus <sub cmd> <args>"

USER root

RUN python3.8 -m pip install --upgrade pip
RUN python3.8 -m pip install numpy==1.21.0 matplotlib==3.3.2 scipy==1.5.2 networkx==2.5.1 biopython==1.76 prody==1.10.11 MDAnalysis==1.1.1 numba==0.53.1
RUN python3.8 -m pip install correlationplus==0.2.0
RUN python3.8 -m pip install ipython

RUN useradd -m correlationplus
USER correlationplus
WORKDIR /home/correlationplus

ENTRYPOINT ["/bin/sh", "-c" , "correlationplus $0 $@"]
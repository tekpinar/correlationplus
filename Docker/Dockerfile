##############################################################################
# correlationplus - A Python package to calculate, visualize and analyze      #
#                   correlation maps of proteins.                             #
# Authors: Mustafa Tekpinar                                                   #
# Copyright (C) Mustafa Tekpinar, 2017-2018                                   #
# Copyright (C) CNRS-UMR3528, 2019                                            #
# Copyright (C) Institut Pasteur Paris, 2020-2021                             #
# Copyright (C) Mustafa Tekpinar, 2022-2026                                   #
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
FROM python:3.11.5


LABEL org.label-schema.vcs-url="https://github.com/tekpinar/correlationplus"
LABEL org.label-schema.version="0.2.3"
LABEL org.label-schema.description="A Python package to calculate, visualize and analyze dynamical correlations of proteins."
LABEL org.label-schema.docker.cmd="docker run -v ~:/home/correlationplus correlationplus <sub cmd> <args>"

USER root

RUN python3.11 -m pip install --upgrade pip
RUN python3.11 -m pip install numpy==1.26.4 matplotlib==3.8.4 scipy==1.9 networkx==2.8 biopython==1.79 prody==2.6.1 MDAnalysis==2.4 numba==0.56
RUN python3.11 -m pip install correlationplus==0.2.3
RUN python3.11 -m pip install ipython

RUN useradd -m correlationplus
USER correlationplus
WORKDIR /home/correlationplus

ENTRYPOINT ["/bin/sh", "-c" , "correlationplus $0 $@"]

# If run command on Mac, you can remove sudo
# Build command
# sudo docker build -t structuraldynamicslab/correlationplus:latest -t structuraldynamicslab/correlationplus:v0.2.3 .
# sudo docker login
# sudo docker push structuraldynamicslab/correlationplus:v0.2.3 

# Pull command on Mac OS
# sudo docker pull structuraldynamicslab/correlationplus:v0.2.3

# Run commands
# sudo docker run -i -t -v $(pwd) --platform linux/amd64 structuraldynamicslab/correlationplus:v0.2.3 /bin/bash

# Bind your current folder to a particular folder in docker.
# In this way, you fall into the bash terminal of the docker image. 
# docker run -ti --rm --mount type=bind,source=$PWD,target=/mnt/data/software structuraldynamicslab/correlationplus:v0.2.3

#If you want to use the docker image like an executable file:
# Please note that you must have the necessary files inside your current folder and run this command!
# sudo docker run --rm -v $PWD:/mnt/data/software structuraldynamicslab/correlationplus:v0.2.3 

FROM python:3.8.5


LABEL org.label-schema.vcs-url="https://github.com/tekpinar/correlationplus"
LABEL org.label-schema.version="0.1.6"
LABEL org.label-schema.description="A Python package to calculate, visualize and analyze dynamical correlations maps of proteins."
LABEL org.label-schema.docker.cmd="docker run -v ~:/home/correlationplus correlationplus <sub cmd> <args>"

USER root

RUN pip3 install numpy==1.19.1 matplotlib==3.3.0 scipy==1.5.2 networkx==2.4 biopython==1.77 prody==1.10.11
RUN pip3 install correlationplus==0.1.6
RUN pip3 install ipython

RUN useradd -m correlationplus
USER correlationplus
WORKDIR /home/correlationplus

ENTRYPOINT ["/bin/sh", "-c" , "correlationplus $0 $@"]
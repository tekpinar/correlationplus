FROM python:3.8.5

LABEL org.label-schema.vcs-url="https://github.com/tekpinar/correlationPlus"
LABEL org.label-schema.version="1.0.3"
LABEL org.label-schema.description="A Python API to plot and analyze dynamical correlations of proteins."
LABEL org.label-schema.docker.cmd="docker run -v ~:/home/correlationPlus correlation_plus <sub cmd> <args>"

USER root

RUN pip3 install numpy==1.19.1 matplotlib==3.3.0 scipy==1.5.2 networkx==2.4 biopython==1.77 prody==1.10.11
RUN pip3 install correlationPlus==0.1.3
RUN useradd -m correlationPlus

USER correlationPlus
WORKDIR /home/correlationPlus

ENTRYPOINT ["/bin/sh", "-c" , "correlationPlus $0 $@"]
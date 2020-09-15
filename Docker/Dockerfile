FROM python:3.8.5

USER root

RUN pip3 install numpy==1.19.1 matplotlib==3.3.0 scipy==1.5.2 networkx==2.4 biopython==1.77 prody==1.10.11
RUN pip3 install -i https://test.pypi.org/simple/ correlationPlus==0.1.3rc2
RUN useradd -m correlationPlus

USER correlationPlus
WORKDIR /home/correlationPlus

ENTRYPOINT ["/bin/sh", "-c" , "correlationPlus $0 $@"]
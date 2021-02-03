Installation
============
Basic Installation
------------------
We recommend installing correlationplus with pip or conda for regular users:

with pip
~~~~~~~~

We recommend to use pip as follows:

``pip install correlationplus``

or if you do not have administration rights

``pip install --user correlationplus``

If you prefer to use a virtualenv
``python3 -m venv correlationplus``
``cd correlationplus``
``source bin/activate``
``pip install correlationplus``


with conda
~~~~~~~~~~

You can also install correlationplus with conda as follows:

``conda install -c bioconda correlationplus``
		
Most of the time, at least one these methods will be sufficient for the installation.
However, if these two methods didn't work for any reason, you can take a look 
to the 'Advanced Installation' instructions.

Advanced Installation
---------------------
If standard installation procedure didn't work for you for any reason, you can 
try one of the methods detailed below:

for developers
~~~~~~~~~~~~~~

We recommend to use pip and a virtualenv
``python3 -m venv correlationplus``
``cd correlationplus``
``source bin/activate``
``mkdir src``
``cd src``
``git clone https://github.com/tekpinar/correlationplus.git # or git@github.com:tekpinar/correlationplus.git``
``cd correlationplus``
``pip install -e .``

from Docker image
~~~~~~~~~~~~~~~~~

Docker images are also available from `Docker Hub <https://hub.docker.com/r/structuraldynamicslab/correlationplus>`_

The computation inside the container will be performed under **correlationplus** id in **/home/correlationplus** directory.
So before running a **correlationplus** container,
do not forget to create and mount a shared directory in the container. 

This directory must be writable.

``mkdir shared_dir``

``cp 6fl9_just_prot_anm_100_modes_rc_15_cross-correlations.txt 6fl9_centeredOrientedAligned2Z.pdb shared_dir``

``chmod 777 shared_dir``

``cd shared_dir``

``docker run -v $PWD:/home/correlationplus structuraldynamicslab/correlation_plus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi``

It is also possible to run an ipython interactive session
``docker run -v $PWD:/home/correlationplus --entrypoint /bin/bash -it structuraldynamicslab/correlationplus:0.1.4rc2``

then once in the container

``ipython``

from Singularity image
~~~~~~~~~~~~~~~~~~~~~~

As the docker image is registered in dockerhub you can also use it directly with `Singularity <https://sylabs.io/docs/>`_

``singularity run docker://structuraldynamicslab/correlationplus diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi``

or in 2 steps

``singularity pull correlationplus.simg docker://structuraldynamicslab/correlation_plus``

``./correlationplus.simg diffMap -i 6fl9_rc15_scalCoeff1_100_modes_lmi_v2.dat -j zacharias_rc15_scalCoeff15_100_modes_lmi.dat -p 6fl9_centeredOrientedAligned2Z.pdb -t lmi``

Unlike docker you do not have to worry about shared directory, your *home* and */tmp* are automatically shared.
You can also run an *ipython* interactive session.
``singularity shell correlationplus.simg``

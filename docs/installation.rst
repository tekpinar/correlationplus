Installation
============

Basic Installation
------------------
We recommend installing correlationplus with pip or conda for regular users:

with pip
~~~~~~~~
We **strongly** encourage you to install correlationplus in a virtualenv::

	python3 -m venv correlationplus
	cd correlationplus
	source bin/activate

Then, you can upgrade pip and install it as follows::

	python3 -m pip install -U pip
	python3 -m pip install correlationplus


or if you do not have administration rights::

	python3 -m pip install --user -U pip
	python3 -m pip install --user correlationplus

If you want to install it without using a virtualenv and encounter an error related to llvmlite, 
you can solve it as follows::

	python3 -m pip install llvmlite --ignore-installed
	python3 -m pip install correlationplus

with conda
~~~~~~~~~~

You can also install correlationplus with conda as follows::

    conda install -c bioconda correlationplus
    
Most of the time, at least one these methods will be sufficient for the installation.
However, if these two methods didn't work for any reason, you can take a look 
at the 'Advanced Installation' instructions to install it from the source.

Advanced Installation
---------------------
If standard installation procedure didn't work for you for any reason, you can 
try one of the methods detailed below:

for developers
~~~~~~~~~~~~~~
The pip version required by some dependencies is >= 21.0.1, which is not the pip version bundled with python 3.(6,7,8)
So, you have to update pip before installing *correlationplus*. Otherwise, you will have trouble during *MDAnalysis* dependency installation.
For this reason, we **strongly** encourage you to install correlationplus in a `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ ::

	python3 -m venv correlationplus
	cd correlationplus
	source bin/activate
	python3 -m pip install -U pip
	python3 -m pip install correlationplus

If you want to install the latest version from the source ::

	python3 -m venv correlationplus
	cd correlationplus
	source bin/activate
	python3 -m pip install -U pip
	mkdir src
	cd src
	git clone https://github.com/tekpinar/correlationplus.git # or git@github.com:tekpinar/correlationplus.git``
	cd correlationplus
	pip install -e .

from Docker image
~~~~~~~~~~~~~~~~~

If you don't want to install or can't install correlatioplus for any reason, you can use docker images. These images do not require any installation. 
You can find docker images for correlationplus  at `<https://hub.docker.com/r/structuraldynamicslab/correlationplus>`_

The computation inside the container will be performed under **correlationplus** id in **/home/correlationplus** directory.
So before running a **correlationplus** container,
do not forget to create and mount a shared directory in the container. 

This directory must be writable to this user. So, you have two possibilities:

1. Map your id on the host to the *correlationplus* user in the container::

       ``-u $(id -u ${USER}):$(id -g ${USER})``
2. Make the shared directory writable to anyone::

       ``chmod 777 shared_dir``

option 1

.. code-block:: shell

    mkdir shared_dir
    cp 6lu7_dimer_no_ligand_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1_ca.pdb shared_dir
    cd shared_dir
    docker run -v $PWD:/home/correlationplus -u $(id -u ${USER}):$(id -g ${USER}) structuraldynamicslab/correlation_plus diffMap\
    				-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
				-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
				-p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

option 2

.. code-block:: shell

    mkdir shared_dir
    cp 6lu7_dimer_no_ligand_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1_ca.pdb shared_dir
    chmod 777 shared_dir
    cd shared_dir
    docker run -v $PWD:/home/correlationplus structuraldynamicslab/correlation_plus diffMap \
    						-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
						-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
						-p 6lu7_dimer_with_N3_protein_sim1-lmi.dat -t lmi


It is also possible to run an ipython interactive session::

    docker run -v $PWD:/home/correlationplus --entrypoint /bin/bash -it structuraldynamicslab/correlationplus:0.1.4rc2

then once in the container

``ipython``

from Singularity image
~~~~~~~~~~~~~~~~~~~~~~

As the docker image is registered in dockerhub you can also use it directly with `Singularity <https://sylabs.io/docs/>`_ ::

    singularity run docker://structuraldynamicslab/correlationplus diffMap \
    					-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
					-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
					-p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

or in 2 steps ::

    singularity pull correlationplus.simg docker://structuraldynamicslab/correlation_plus
    ./correlationplus.simg diffMap \
    			-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
			-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
			-p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

Unlike docker, you do not have to worry about shared directory, your *home* and */tmp* are automatically shared.
You can also run an *ipython* interactive session ::

    singularity shell correlationplus.simg

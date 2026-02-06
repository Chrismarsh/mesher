Development env
-------------------

.. warning::
   Using mesher this way (i.e., built from source) will require setting the ``mesher`` binary path as described in :ref:`Configuration:Environment Variables`.


Building mesher from source has two parts 

1) compile the backend C++ binary and
2) setup a python environment for the python frontend.

Building the dev environment via spack is the supported method.

Steps
***********

An out of source build is recommended.

::

    cd ~/
    git clone https://github.com/Chrismarsh/mesher.git
    cd mesher
    spack env activate .
    spack install # builds dependencies


The backend mesher binary can be independently built by calling cmake in an out of source build.

::

    mkdir build-backend
    cd build-backend
    cmake ../
    make

Usage
*******

::

 PYTHONPATH=~/code/mesher/;MESHER_EXE=~/mesher/build-backend/build-bin/mesher python -m src.mesher.mesher_cli granger_veg.py



Deployment
==========

Use the spack environment defined in ``publish.sh``









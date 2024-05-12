Compilation
-------------

.. warning::
   Using mesher this way (i.e., built from source) will require setting the ``mesher`` binary path as described in :ref:`Configuration:Environment Variables`.


Building mesher from source has two parts 

1) compile the backend C++ binary and
2) setup a python environment for the python frontend.

Building the dev environment via spack is the supported method.

Clone repo
***********

An out of source build is recommended.

::

    cd ~/
    git clone https://github.com/Chrismarsh/mesher.
    cd mesher
    spack env create . py-mesher-spack.yaml
    spack env activate .
    spack install # builds dependencies


The backend mesher binary can be independently built by calling cmake in an out of source build.


What to do if things aren't working
=====================================

On Macos, homebrew tends to update often. Thus, python packages that are installed against homebrew libraries may break in interesting ways.
The most common issue is that the pygdal package hasn't been update for the newest version of gdal in homebrew. This will manifest as an error like
``gdal 2.3.4 != gdal 3.2.1`` . Unfortunately there is not much that can be until pygdal gets updated. It is best to open an issue either on the mesher or nextgis/pygdal githubs.

If there is an install error about ``gdal-config`` missing, please ensure that gdal is installed from your package manager and is available on the path. Running
``gdal-config --version`` should produce output


Deployment
==========

Use the spack environment defined in ``publish.sh``









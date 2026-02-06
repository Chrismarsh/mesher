Installation
============

.. :: warning
    Mesher is only supported on Macos and Linux

Installation of mesher is possible via ``pip``, ``conda``, or ``spack``. A ``spack``  build will be the easiest due to
all the potentially complex dependencies of vtk and gal being handled automatically. See below on how to install with
``pip`` which requires system library installs. ``conda`` builds are particularly tested and YMMV.


Wheels are not prebuilt for mesher. Instead, mesher will need to be compiled as part of the install step.
This thus requires a functional build environment.


spack
+++++++

 ``py-mesher`` has not been added to spack's built in repos yet, so please add it's spack repo.

- ``spack repo add https://github.com/Chrismarsh/spack-repo.git``
- ``spack install py-mesher``


conda
++++++

::
    # ensure conda-forge is added
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    conda install mesher


pip
++++

Using pip, mesher must be built against system libraries.

.. note::
   Depending on your python install, ``pip`` may be ``pip3``


Ensure the following are installed via package manager:

For macos:

::

      brew install gdal
      brew install boost
      brew install cgal
      brew install metis

For Ubuntu:

::

    sudo apt-get install libgdal-dev
    sudo apt-get install libcgal-dev
    sudo apt-get install gdal-bin
    sudo apt-get install libcgal-dev
    sudo apt-get install libboost-filesystem-dev
    sudo apt-get install libboost-program-options-dev
    sudo apt-get install libmetis-dev

    # on Ubuntu 20.04+
    sudo apt-get install python3-gdal
    # prior to Ubuntu 20.04, use this instead of python3-gdal
    # sudo apt-get install python-gdal


.. :: warning
    On linux you may need ``libffi`` if, upon running ``pip``, there is an error about ``_ctypes``

    On Ubuntu
    ``apt-get install libffi-dev``

    On CentOS/Fedora
     ``dnf install libffi-devel``

Then install mesher with

::

    pip install mesher


















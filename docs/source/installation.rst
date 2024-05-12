Installation
============

.. :: warning
    Mesher is only supported on Macos and Linux

Installation of mesher is possible via ``pip``, ``conda``, or ``spack``. A ``conda`` build will be the easiest due to
all the potentially complex dependencies of vtk and gal being handled automatically. See below on how to install with
``pip`` which requires system library installs. In a HPC or development context, ``spack`` is prefered.


Wheels are not prebuilt for mesher. Instead, mesher will need to be compiled as part of the install step.
This thus requires a functional build environment.

It is easiest if a Python version that has vtk wheels is used -- see `here <https://pypi.org/project/vtk/9.1.0/#files:vtk>`_ for details on wheel availability.
Consider using the `pyenv <https://github.com/pyenv/pyenv>`_ python version manager if this is not your system default Python version.

conda
++++++

::
    # ensure conda-forge is added
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    conda install mesher

spack
+++++++

 ``py-mesher`` has not been added to spack's built in repos yet, so it needs to be added.

- Clone https://github.com/Chrismarsh/spack-repo
- Add `spack-repo` to spack `repos.yaml` https://spack.readthedocs.io/en/latest/repositories.html
- `spack install py-mesher`


pip
++++

Using pip, mesher must be built against system libraries.

.. note::
   Depending on your python install, ``pip`` may be ``pip3``

System
--------

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




Install of github branch
++++++++++++++++++++++++++
You can optionally use pip to install the most recent github version or a github branch. However, the automatic
setup of the build environment does not occur, so ensure ``scikit-build-core`` is installed. Then,

::

    pip install git+https://github.com/Chrismarsh/mesher@branch-name

If mesher is already installed, use ``--force-reinstall`` to reinstall it.

















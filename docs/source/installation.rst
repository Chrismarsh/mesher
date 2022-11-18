Installation
============

Installation of mesher is possible via ``pip``.  Installation into a conda environment probably works but is not tested.
Mesher is only supported on Macos and Linux for Python 3.7+

Wheels are not prebuilt for mesher. Instead, mesher will need to be compiled as part of the pip install step. This thus requires a functional build environment.

Mesher is tested on Macos (brew) and Ubuntu (apt-get) although other configurations likely work as expected. Adjust the dependencies below as needed.

It is easiest if Python 3.7, 3.8, 3.9 is used as one of the dependencies (vtk) has prebuilt wheels (see `here <https://pypi.org/project/vtk/9.1.0/#files:vtk>`_ for details on wheel availability).
Consider using the `pyenv <https://github.com/pyenv/pyenv>`_ python version manager if this is not your system default Python version.

Setup environment
+++++++++++++++++++

Mesher can be built against system libraries or against conan libraries.

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

Conan
--------
Install conan via

::

    pip install conan

And then setup a new profile

::

    conan profile new default --detect
    conan config install https://github.com/Chrismarsh/conan-config.git


This configuration file setups use of revisions, two new remotes (bincrafters, CHM), and tweaks the ``settings.yml`` file to have ubuntu-18.04 and ubuntu-20.04 distros. Setting
``os.distro = 'ubuntu-20.04'`` will enable the use of prebuilt library binaries.

Then setup conan to use the new C++ ABI and C++ standard

::

  conan profile update settings.compiler.cppstd=14 default

If using clang (e.g.,Macos), do

::

   conan profile update settings.compiler.libcxx=libc++ default  #with clang

and if using gcc, do

::

   conan profile update settings.compiler.libcxx=libstdc++11 default  #with gcc

then install mesher with

::

    USE_CONAN=TRUE pip install mesher


conda
++++++

.. warning::
    This is not tested! Mixing conan + conda seems to not be reliable so please use system libraries.

The Anaconda python environment supports ``pip`` installs. This example shows installing Anaconda, however if you already have Anaconda installed, then only the instructions from ``conda create`` onward is required.

::

  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/conda
  source $HOME/conda/bin/activate
  conda init
  conda update -y --all
  conda create -y --name mesher python=3.7
  conda activate mesher
  pip install mesher

This approach will use the system installed gdal.



Install of github branch
++++++++++++++++++++++++++
You can optionally use pip to install the most recent github version or a github branch. However, the automatic
setup of the build environment does not occur, so ensure ``scikit-build``, ``cmake``, ``conan``, and ``ninja`` are installed. Then,

::

    pip install git+https://github.com/Chrismarsh/mesher@branch-name

If mesher is already installed, use ``--force-reinstall`` to reinstall it.

















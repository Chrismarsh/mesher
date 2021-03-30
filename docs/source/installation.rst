Installation
============

Installation of mesher is possible via ``pip``. 

Mesher is only supported on Macos and Linux for Python 3.7 and 3.8

::

   pip install mesher


It is easiest if Python 3.7, or 3.8 is used (see :ref:`build:vtk` for details on wheel availability).
Consider using the `pyenv <https://github.com/pyenv/pyenv>`_ python version manager as described in the :ref:`build:Setup Python` section. 


Ensure the following are installed via package manager:

For macos:

::

      brew install gdal
      brew install boost
      brew install cgal

For Ubuntu:

::

    sudo apt-get install libgdal-dev
    sudo apt-get install python-gdal
    sudo apt-get install gdal-bin
    sudo apt-get install libcgal-dev
    sudo apt-get install libboost-filesystem-dev
    sudo apt-get install libboost-program-options-dev

On Ubuntu 20.04, use 

::

   sudo apt-get install python3-gdal

.. :: warning
    On linux you may need ``libffi`` if, upon running ``pip``, there is an error about ``_ctypes``

    On Ubuntu
    ``apt-get install libffi-dev``

    On CentOS/Fedora
     ``dnf install libffi-devel``

.. note::
   If ``conan`` is used during the compilation, two new remotes will be automatically added to the ``.conan/remotes.json`` file:
   ::

      conan remote add bincrafters https://api.bintray.com/conan/bincrafters/public-conan
      conan remote add CHM https://api.bintray.com/conan/chrismarsh/CHM

   These are needed to download the required dependencies for the backend

.. note::
   Depending on your python install, ``pip`` may be ``pip3``

Full working example
**********************

This provides a full example of setting up a python environment with the assumption the user does not have a python environment setup.

.. note::
   If on macos, and you have homebrew, you can optionally install ``pyenv`` and ``pyenv-virtualenv`` via brew.


Install `pyenv`_
::

   curl https://pyenv.run | bash

Note, this automatically installs ``pyenv-virtualenv``.

Restart the shell and update
::

   exec $SHELL 
   pyenv update


Install python 3.7.6
::

   pyenv install 3.7.6
   pyenv shell 3.7.6

Create new virtual environment for ``mesher`` and install mesher
::
   
   pyenv virtualenv 3.7.6 mesher-3.7.6
   pyenv activate mesher-3.7.6
   pip install mesher


Now, if you wish to run mesher activate the virtualenv
::
   
   pyenv activate mesher-3.7.6
   mesher.py my-config-file.py


conda
******
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
*************************
You can optionally use pip to install the most recent github version or a github branch. However, the automatic
setup of the build environment does not occur, so ensure ``scikit-build``, ``cmake``, ``conan``, and ``ninja`` are installed. Then,

::

    pip install git+https://github.com/Chrismarsh/mesher@branch-name

If mesher is already installed, use ``--force-reinstall`` to reinstall it.


Automatic virtualenv activation
*******************************

The automatic virtualenv activation provided by ``pyenv-virtualenv`` can make it easier to work with virtual environments. 

Follow point 2 `here <https://github.com/pyenv/pyenv-virtualenv>`_ to enable this feature.

Any folder with a ``.python-version`` that contains a  valid virtualenv specification will have it automatically enabled upon entering that folder. For example,

::
   
   cd my-working-folder
   echo "mesher-3.7.6" >> .python-version


will automatically activate the above-created virtualenv every time that folder is entered, and deactivate when leaving.


















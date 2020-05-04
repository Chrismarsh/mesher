Installation
============

Installation of mesher is possible via ``pip``. 

Mesher is only supported on Macos and Linux for Python >= 3.6

::

   pip install mesher


It is easiest if Python 3.6, or 3.7 is used (see :ref:`build:vtk` for details on wheel availability). 
Consider using the `pyenv <https://github.com/pyenv/pyenv>`_ python version manager as described in the :ref:`build:Setup Python` section. 


.. warning::
   Ensure the following are installed via package manager:
      - gdal development libaries (e.g., ``dnf install gdal-devel`` or ``brew install gdal``). See other notes :ref:`build:gdal` as you may also need the binary package.
      - vtk development libaries (e.g., ``dnf install vtk-devel`` or ``brew install vtk``) **if not** Python 3.6, or 3.7



.. note::
   As ``conan`` is used during the compilation, two new remotes will be automatically added to the ``.conan/remotes.json`` file:
   ::

      conan remote add bincrafters https://api.bintray.com/conan/bincrafters/public-conan
      conan remote add CHM https://api.bintray.com/conan/chrismarsh/CHM

   These are needed to download the required dependencies for the backend

.. note::
   Depending on your python install, ``pip`` may be ``pip3``

Full working example
**********************

This provides a full example of setting up a python environemtn with the assumption the user does not have a python environment setup. 

.. note::
   If on macos, and you have homebrew, install ``pyenv`` and ``pyenv-virtualenv`` via brew. 


Install `pyenv`_
::

   curl https://pyenv.run | bash

Install `pyenv-virtualenv <https://github.com/pyenv/pyenv-virtualenv>`_ 
::
   
   git clone https://github.com/pyenv/pyenv-virtualenv.git $(pyenv root)/plugins/pyenv-virtualenv

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





Automatic virtualenv activation
*******************************

The automatic virtualenv activation provided by ``pyenv-virtualenv`` can make it easier to work with virtual environments. 

Follow point 2 `here <https://github.com/pyenv/pyenv-virtualenv>`_ to enable this feature.

Any folder with a ``.python-version`` that contains a  valid virtualenv specification will have it automatically enabled upon entering that folder. For example,

::
   
   cd my-working-folder
   echo "mesher-3.7.6" >> .python-version


will automatically activate the above-created virtualenv every time that folder is entered, and deactivate when leaving.


















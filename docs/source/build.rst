Compilation
-----------

.. warning::
   Using mesher this way (i.e., built from source) will require setting the ``mesher`` binary path as described in :ref:`Configuration:Environment Variables`.


Build from source
=================

Building mesher from source has two parts 

1) compile the backend C++ binary and
2) setup a python environment for the python frontend.


Build Requirements
*******************

Linux and MacOS are the only supported environments.
Python  3.7 or 3.8

Python 3.9 does not have vtk wheels

- cmake >= 3.16
- C++11 compiler (gcc 7.x+ recommended)

On MacOS, homebrew should be used to install cmake and conan. Macport based installs likely work, but have not been tested.

Compile
********

Building mesher requires:

- cmake (>3.16)
- boost (>1.70)
- gdal (>2.4 or >3.0)
- cgal


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

Throughout, this document assumes a working development environment

Clone repo
***********
An out of source build is recommended. For example:
::

    cd ~/
    git clone https://github.com/Chrismarsh/mesher.
    mkdir ~/build && cd ~/build


Run cmake
*********
You can set the install prefix to be anywhere, such as shown in the example below
::

    cmake ~/mesher -DCMAKE_INSTALL_PREFIX=/opt/mesher

This should complete without any errors.


Building
*********
Using make
::

    make install 


Setup Python
============


The python vtk bindings (see below) are most easily installed using the `wheels packages <https://pypi.org/project/vtk/#files>`_. However, vtk wheels only exist for Python 3.6, and 3.7.
Therefore it's highly recommended to use Python 3.7. Doing so can easily be done with `pyenv <https://github.com/pyenv/pyenv>`_ to manage python versions:
::

   pyenv install 3.7.6
   pyenv shell 3.7.6 # activate this version of python for this shell


If ``pyenv`` is used, then the excellent `pyenv-virtualenv <https://github.com/pyenv/pyenv-virtualenv>`_ wrapper can easily streamline ``virtualenv`` creation 
::

   pyenv virtualenv 3.7.6 mesher-3.7.6
   pyenv activate mesher-3.7.6


Regardless of how the virtualenv is setup, install the following:


vtk
***

vtk `wheels <https://prabhuramachandran.blogspot.com/2018/01/vtk-810-wheels-for-all-platforms-on-pypi.html>`_ only `exist <https://pypi.org/project/vtk/#files>`_ for Python 3.5, 3.6, and 3.7. If building from source, ensure vtk development files (e.g., `vtk-devel`) are installed through your system's package manager.

::

   pip install vtk


gdal 
****

It's recommended that gdal python bindings are installed via `pygdal <https://github.com/nextgis/pygdal>`_. gdal doesn't provide wheels, so ``pygdal`` will need to build from source. Therefore ensure gdal development files (e.g., ``gdal-devel``) are installed through your system's package manager. 

On linux, depending on the distro used, you may need to also install the gdal binaries. On Ubuntu this is
::

   sudo apt-get install libgdal-dev
   sudo apt-get install gdal-bin

The system gdal binaries such as ``gdalwarp`` are only avilable in the ``gdal-bin`` pacakage.


other libraries
***************
:: 
   
   pip install numpy scipy cloudpickle



Deployment
==========
Notes for how to deploy to Pypi:

:: 
   
   pip install scikit-build
   pip install twine
   pip install wheel

::

   python setup.py sdist bdist_wheel
   twine upload  dist/*


Note that version number needs to be incremented for each Pypi upload


Optional Conan
===============
Optionall mesher's dependencies can be build using `conan <https://www.conan.io/>`_ for dependency management.

All of the mesher dependencies are built on Github-CI and uploaded to the bintray repository to serve prebuilt binaries. This means that if the mesher build is done with supported compilers and operating system (described later), the dependencies do not need to be built by the end user.


.. warning::
   The python gdal bindings uses a system-wide gdal rather than the conan gdal the mesher C++ backend links against. This will hopefully be resolved in the future. However, as no data passes between the C++ and Python, having different gdal versions poses no problem.


Setup dependencies
******************
Install the dependencies into your local conan cache (`~/.conan/data`)
::

    cd ~/build #if you have not already
    conan install ~/mesher -if=.


The `-if=.` will produce the ``FindXXX.cmake`` files required for the mesher build in the current directory.

If you need to build dependencies from source, use the `--build missing` option like:
::

    conan install ~/mesher -if=. --build missing

Then enable the conan build in the CMake file

::

    cmake ~/mesher -DCMAKE_INSTALL_PREFIX=/opt/mesher -DUSE_CONAN=True














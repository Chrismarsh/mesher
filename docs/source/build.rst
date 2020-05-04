Compilation
-----------

Build from source
=================

Building mesher from source has two parts 

1) compile the backend C++ binary and
2) setup a python environment for the python frontend.

The easiest way to build mesher is to use `conan <https://www.conan.io/>`_ for dependency management. 

All of the mesher dependencies are built on Travis-CI and uploaded to the bintray repository to serve prebuilt binaries. This means that if the mesher build is done with supported compilers and operating system (described later), the dependencies do not need to be built by the end user.

.. warning::
   Using mesher this way will require setting the ``mesher`` binary path as described in :ref:`Configuration:Environment Variables`.


Build Requirements
*******************

Linux and MacOS are the only supported environments.
Python >= 3.6

If bintray binaries can be used, the only requirements are:

- conan >= 1.21
- cmake >= 3.16
- C++11 compiler (gcc 7.x+ recommended)

On MacOS, homebrew should be used to install cmake and conan. Macport based installs likely work, but have not been tested.

Compile
********

Building mesher requires:

- cmake (>3.16)
- boost (>1.70)
- cgal


Throughout, this document assumes a working development environment, but a blank conan environment.

Setup conan
***********

Initialize a new profile
::

    conan profile new default --detect
    conan remote add bincrafters https://api.bintray.com/conan/bincrafters/public-conan
    conan remote add CHM https://api.bintray.com/conan/chrismarsh/CHM
    conan profile update settings.compiler.cppstd=14 default  


conan needs to be told to use new C++11 ABI. If using clang (e.g., MacOs), do
::

    conan profile update settings.compiler.libcxx=libc++ default  #with clang


and if using gcc, do
::

    conan profile update settings.compiler.libcxx=libstdc++11 default  #with gcc


If you change compilers, such as on a cluster with a modules system, you can rerun 
::
    
    conan profile new default --detect --force


to detect the new compiler settings. The ``cppstd`` and ``libcxx`` settings need to be reapplied once this is done.

Intel compiler
**************
If you're using the Intel compiler, ensure ``compilervars.sh`` is sourced, e.g.,
::

    source /opt/intel/bin/compilervars.sh intel64

prior to running conan. Use the above gcc settings for conan.


Clone repo
***********
An out of source build is recommended. For example:
::

    cd ~/
    git clone https://github.com/Chrismarsh/mesher.
    mkdir ~/build && cd ~/build

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

Run cmake
*********
You can set the install prefix to be anywhere, such as shown in the example below
::

    cmake ~/mesher -DCMAKE_INSTALL_PREFIX=/opt/mesher


This should complete without any errors.


If using the Intel compiler, add the following cmake flags:
::

    -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_FORTRAN_COMPILER=ifort

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

.. note::
   The python gdal bindings uses a system-wide gdal rather than the conan gdal the mesher C++ backend links against. This will hopefully be resolved in the future. However, as no data passes between the C++ and Python, having different gdal versions poses no problem.

On linux, depending on the distro used, you may need to also install the gdal binaries and, paradoxically, the gdal python bindings. On Ubuntu this is
::

   sudo apt-get install libgdal-dev
   sudo apt-get install gdal-bin
   sudo apt-get install python-gdal

The system gdal-python bindings is because certain python scripts such as ``gdal_polygonize.py`` are only available when installing ``python-gdal``, and the gdal binaries such as ``gdalwarp`` are only avilable in the ``gdal-bin`` pacakage.


other libraries
***************
:: 
   
   pip install numpy scipy



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




















Compilation
-------------

.. warning::
   Using mesher this way (i.e., built from source) will require setting the ``mesher`` binary path as described in :ref:`Configuration:Environment Variables`.


Building mesher from source has two parts 

1) compile the backend C++ binary and
2) setup a python environment for the python frontend.


Setup the build environment as described on the install page

Clone repo
***********

An out of source build is recommended. For example:

::

    cd ~/
    git clone https://github.com/Chrismarsh/mesher.
    mkdir ~/build && cd ~/build


Python packages
********************

::

    setuptools
    wheel
    conan
    scikit-build>=0.11.1
    ninja
    vtk
    pygdal-chm=="`gdal-config --version`.*"
    numpy
    scipy
    matplotlib
    cloudpickle
    metis

(Optional) Conan dependencies
*********************************
Optionally mesher's dependencies can be build using `conan <https://www.conan.io/>`_ for dependency management.

All of the mesher dependencies are built on Github-CI and uploaded to the bintray repository to serve prebuilt binaries. This means that if the mesher build is done with supported compilers and operating system (described later), the dependencies do not need to be built by the end user.

.. warning::
   The python gdal bindings uses a system-wide gdal rather than the conan gdal the mesher C++ backend links against. This will hopefully be resolved in the future. However, as no data passes between the C++ and Python, having different gdal versions poses no problem.

.. warning::
    Conan and conda don't seem to consistently work. Use at your own risk.


Setup Conan as described on the installation page.

Install the dependencies into your local conan cache (`~/.conan/data`)

::

    cd ~/build #if you have not already
    conan install ~/mesher -if=. --build missing


The `-if=.` will produce the ``FindXXX.cmake`` files required for the mesher build in the current directory, building missing as needed.


Build
***********

You can set the install prefix to be anywhere, such as shown in the example below

::

    cmake ~/mesher -DCMAKE_INSTALL_PREFIX=/opt/mesher

if you're using conan,

::

    cmake ~/mesher -DCMAKE_INSTALL_PREFIX=/opt/mesher -DUSE_CONAN=True

then

::

    make install


What to do if things aren't working
=====================================

On Macos, homebrew tends to update often. Thus, python packages that are installed against homebrew libraries may break in interesting ways.
The most common issue is that the pygdal package hasn't been update for the newest version of gdal in homebrew. This will manifest as an error like
``gdal 2.3.4 != gdal 3.2.1`` . Unfortunately there is not much that can be until pygdal gets updated. It is best to open an issue either on the mesher or nextgis/pygdal githubs.

If there is an install error about ``gdal-config`` missing, please ensure that gdal is installed from your package manager and is available on the path. Running
``gdal-config --version`` should produce output


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









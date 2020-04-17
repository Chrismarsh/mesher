Compilation
-----------

Background
=============
The easiest way to build mesher is to use `conan <https://www.conan.io/>`_ for dependency management. 

All of the mesher dependencies are built on Travis-CI and uploaded to the bintray repository to serve prebuilt binaries. This means that if the mesher build is done with supported compilers and operating system (described later), the dependencies do not need to be built by the end user.

Build Requirements
===================
Linux and MacOS are the only supported environments.

If bintray binaries can be used, the only requirements are:

- conan >= 1.21
- cmake >= 3.16
- C++11 compiler (gcc 7.x+ recommended)

On MacOS, homebrew should be used to install cmake and conan. Macport based installs likely work, but have not been tested.

Compile
========
Building mesher requires:

- cmake (>3.16)
- boost (>1.70)
- cgal

Throughout, this document assumes a working development environment, but a blank conan environment.

Setup conan
===========
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


to detect the new compiler settings. The `cppstd` and `libcxx` settings need to be reapplied once this is done.

Intel compiler
==============
Ensure the Intel compilervars is sourced, e.g.,
::

    source /opt/intel/bin/compilervars.sh intel64

prior to running the conan. Use the above gcc settings for conan.


Clone repo
===========
An out of source build is recommended. For example:
::

    cd ~/
    git clone https://github.com/Chrismarsh/mesher.
    mkdir ~/build && cd ~/build

Setup dependencies
===================
Install the dependencies into your local conan cache (`~/.conan/data`) 
::
    
    cd ~/build #if you have not already
    conan install ~/mesher -if=.


The `-if=.` will produce the `FindXXX.cmake` files required for the mesher build in the current directory. 

If you need to build dependencies from source, use the `--build missing` option like:
::

    conan install ~/mesher -if=. --build missing

Run cmake
=========
You can set the install prefix to be anywhere, such as shown in the example below
::

    cmake ~/mesher -DCMAKE_INSTALL_PREFIX=/opt/mesher


This should complete without any errors.


If using the Intel compiler, add the following cmake flags:
::

    -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_FORTRAN_COMPILER=ifort

Building
========
Using make
::

    make install 































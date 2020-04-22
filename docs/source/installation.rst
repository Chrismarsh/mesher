Installation
============

Installation of mesher is possible via `pip`.

It is easiest if Python 3.5., 3.6, or 3.7 is used (see `build`_ for details on `vtk` wheel availability). 
Consider using the `pyenv <https://github.com/pyenv/pyenv>`_ python version manager as described in the `build`_ section. 

.. note::
   Ensure the following are installed via package manager:
      - gdal development libaries (e.g., `gdal-devel`)
      - vtk development libaries (e.g., `vtk-devel`) **if** not Python 3.5., 3.6, or 3.7

::

   pip install mesher


.. note::
   As `conan` is used during the compilation, two new remotes will be automatically added to the `.conan/remotes.json` file:
   ::

      conan remote add bincrafters https://api.bintray.com/conan/bincrafters/public-conan
      conan remote add CHM https://api.bintray.com/conan/chrismarsh/CHM

   These are needed to download the required dependencies for the backend


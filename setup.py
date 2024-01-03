from skbuild import setup
import subprocess
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version
import packaging.version
from setuptools import find_packages
import os

def gdal_dependency():
    try:
        version = subprocess.run(["gdal-config", "--version"], stdout=subprocess.PIPE).stdout.decode()
        version = version.replace('\n', '')

        gdal_depends = ''
        if packaging.version.parse(version) >= packaging.version.parse("3.5.0"):
            # >= 3.5 required for this type of gdal python binding install
            gdal_depends = f'GDAL[numpy]=={version}.*'
        else:
            gdal_depends = f'pygdal=={version}.*'

        return gdal_depends
    except FileNotFoundError as e:
        raise(""" ERROR: Could not find the system install of GDAL. 
                  Please install it via your package manage of choice.
                """)

# Add CMake as a build requirement if cmake is not installed or is too low a version
# https://scikit-build.readthedocs.io/en/latest/usage.html#adding-cmake-as-building-requirement-only-if-not-installed-or-too-low-a-version
setup_requires = []
try:
    if packaging.version.parse(get_cmake_version()) < packaging.version.parse("3.16"):
        setup_requires.append('cmake')
except SKBuildError:
    setup_requires.append('cmake')


USE_CONAN = False
try:
  USE_CONAN = os.environ["USE_CONAN"]
except KeyError as e:
  pass # it's ok we don't have this

USE_CONAN = str(USE_CONAN).upper() 

setup(name='mesher',
    version='2.0.7',
    description='Landsurface model mesh generation',
    long_description="""
    Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. 
    Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation.
    The paper describing mesher can be `found here <https://www.usask.ca/hydrology/papers/Marsh,_et_al_2018.pdf>`_
    
    Key points
    
    *  A novel multi-objective unstructured mesh generation software, Mesher
    *  Heterogeneity in topography is resolved as well as hydrologically important surface and sub-surface attributes
    *  Spatial heterogeneity is better preserved compared to existing mesh generators
    """,
    author='Chris Marsh',
    author_email='chris.marsh@usask.ca',
    url="https://github.com/Chrismarsh/mesher",
    packages=find_packages(where="pysrc"),
    package_dir={
        '': 'pysrc',
    },
    include_package_data=True,
    cmake_args=['-DCMAKE_BUILD_TYPE:STRING=Release',
              '-DUSE_CONAN:BOOL='+USE_CONAN],
    scripts=["mesher.py","tools/mesh2vtu.py", "tools/meshmerge.py","tools/meshpermutation.py","tools/meshstats.py", "tools/mesh2shp.py"],
    install_requires=['vtk','numpy', 'scipy', 'matplotlib', 'cloudpickle', gdal_dependency(),
                      'metis', 'mpi4py', 'natsort'],
    setup_requires=setup_requires,
    python_requires='>=3.7'
    )
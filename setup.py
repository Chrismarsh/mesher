from skbuild import setup
import subprocess
from packaging.version import LegacyVersion
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version
import os

def get_installed_gdal_version():
    try:
        version = subprocess.run(["gdal-config","--version"], stdout=subprocess.PIPE).stdout.decode()

        version = version.replace('\n', '')
        version = "=="+version+".*"
        return version
    except FileNotFoundError as e:
        raise(""" ERROR: Could not find the system install of GDAL. 
                  Please install it via your package manage of choice.
                """
            )

# Add CMake as a build requirement if cmake is not installed or is too low a version
# https://scikit-build.readthedocs.io/en/latest/usage.html#adding-cmake-as-building-requirement-only-if-not-installed-or-too-low-a-version
setup_requires = []
try:
    if LegacyVersion(get_cmake_version()) < LegacyVersion("3.16"):
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
      version='1.5.6',
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
      include_package_data=True,
      cmake_args=['-DCMAKE_BUILD_TYPE:STRING=Release',
                  '-DUSE_CONAN:BOOL='+USE_CONAN],
      scripts=["mesher.py","tools/mesh2vtu.py", "tools/meshmerge.py","tools/meshpermutation.py","tools/meshstats.py"],
      install_requires=['vtk','pygdal'+get_installed_gdal_version(),'numpy','scipy','matplotlib','cloudpickle'],
      setup_requires=setup_requires,
      python_requires='>=3.7, <3.9'
     )
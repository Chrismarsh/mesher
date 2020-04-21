from skbuild import setup
import os
import subprocess

def get_installed_gdal_version():
    try:
        version = subprocess.check_output(["gdal-config","--version"],text=True)
        version = version.replace('\n', '')
        version = "=="+version+".*"
        return version
    except FileNotFoundError as e:
        raise(""" ERROR: Could not find the system install of GDAL. 
                  Please install it via your package manage of choice.
                """
            )




setup(name='mesher',
      version='1.0',
      description='Landsurface scheme mesh generation',
      author='Chris Marsh',
      author_email='chris.marsh@usask.ca',
      url='https://github.com/Chrismarsh/mesher',

      # cmake_install_dir='./bin ',
      cmake_args=['-DCMAKE_BUILD_TYPE=Release'],
      scripts=["mesher.py","permutation_tools.py"],
      install_requires=['vtk','pygdal'+get_installed_gdal_version(),'numpy','scipy']
     )
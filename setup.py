from skbuild import setup
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
      version='1.0.1',
      description='Landsurface scheme mesh generation',
      long_description="""
      Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. 
      Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation.
      The paper describing mesher can be `found here <https://www.usask.ca/hydrology/papers/Marsh,_et_al_2018.pdf>`_

        Key points
        
        *  A novel multi-objective unstructured mesh generation software, Mesher, is presented
        *  Heterogeneity in topography is resolved as well as hydrologically important surface and sub-surface attributes
        *  Spatial heterogeneity is better preserved compared to existing mesh generators
      """,
      author='Chris Marsh',
      author_email='chris.marsh@usask.ca',
      url="https://github.com/Chrismarsh/mesher",
      include_package_data=True,
      cmake_args=['-DCMAKE_BUILD_TYPE=Release'],
      scripts=["mesher.py","tools/mesh2vtu.py", "tools/meshmerge.py","tools/meshpermutation.py","tools/meshstats.py"],
      install_requires=['vtk','pygdal'+get_installed_gdal_version(),'numpy','scipy']
     )
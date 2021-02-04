# mesher

Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation. The paper describing *mesher* can be [found here](https://www.usask.ca/hydrology/papers/Marsh,_et_al_2018.pdf).

![](docs/source/images/mesher_veg.png)

### How to use
Detailed documentation is given [here](https://mesher-hydro.readthedocs.io).

### Install

Detailed documentation on how to install is given [here](https://mesher-hydro.readthedocs.io/en/latest/installation.html).

In short:
  - Python 3.7, 3.8
  - Install GDAL development files from package manager (e.g., `gdal-devel`)
  - C++14 compliant gcc (gcc 7.3+ is good)
  - You **might** need to install `vtk-devel` from package manager if the python wheels for vtk don't exist on linux.
      + Python 3.9 doesn't have VTK wheels yet on macos
  
Then:

```
$ pip install mesher
```






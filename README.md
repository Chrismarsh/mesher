# mesher

Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation. The paper describing *mesher* can be [found here](https://www.usask.ca/hydrology/papers/Marsh,_et_al_2018.pdf).

![](docs/source/images/mesher_veg.png)

### How to use
Detailed documentation is given [here](https://mesher-hydro.readthedocs.io).

### Install

Detailed documentation on how to install is given [here](https://mesher-hydro.readthedocs.io/en/latest/installation.html).

In short:
  - Python >=3.6
  - Install GDAL development files from package manager (e.g., `gdal-devel`)
  - Have Python 3.x and a C++14 compliant gcc (gcc 7.3+ is good)
  - If you **aren't** using Python 3.6, or 3.7, install `vtk-devel` from package manager
  
Then:

```
$ pip install mesher
```






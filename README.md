# mesher

Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation. The paper describing it can be [found here](https://www.usask.ca/hydrology/papers/Marsh,_et_al_2018.pdf).

### Key points
*	A novel multi-objective unstructured mesh generation software, Mesher, is presented
*	Heterogeneity in topography is resolved as well as hydrologically important surface and sub-surface attributes
*	Spatial heterogeneity is better preserved compared to existing mesh generators

### Documentation

Detailed documentation on how to install is given in the [documentation](https://mesher-hydro.readthedocs.io).

In short:
  - Install GDAL development files from package manager (e.g., `gdal-devel`)
  - Have Python 3.x and a C++14 compliant gcc (gcc 7.3+ is good)
  - If you **aren't** using Python 3.5, 3.6, or 3.7, install `vtk-devel` from package manager
  
Then:

```
$ pip install mesher
```



![](docs/source/images/mesher_veg.png)


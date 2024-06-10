# mesher

Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation. The paper describing *mesher* can be [found here](https://research-groups.usask.ca/hydrology/documents/pubs/papers/marsh,_et_al_2018.pdf).

![](docs/source/images/mesher_veg.png)

### How to use
Detailed documentation is given [here](https://mesher-hydro.readthedocs.io).

### Install
Build requirements
  - Python >= 3.7
  - C++14 compliant gcc (>= gcc 7.3)
  - gdal >=3.8, cgal, boost, vtk>=9, metis

```
$ pip install mesher
```

or

```
$ conda install mesher
```

Detailed documentation on how to install is given [here](https://mesher-hydro.readthedocs.io/en/latest/installation.html).




#### Spack
- Clone https://github.com/Chrismarsh/spack-repo
- Add `spack-repo` to spack `repos.yaml` https://spack.readthedocs.io/en/latest/repositories.html
- `spack install py-mesher`





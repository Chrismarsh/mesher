# mesher

Mesher is a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation. The paper describing it can be [found here](https://www.usask.ca/hydrology/papers/Marsh,_et_al_2018.pdf).

### Key points
*	A novel multi-objective unstructured mesh generation software, Mesher, is presented
*	Heterogeneity in topography is resolved as well as hydrologically important surface and sub-surface attributes
*	Spatial heterogeneity is better preserved compared to existing mesh generators


An example of a mesh that has been generated to represent topography and vegetation is shown below for the Rocky Mountains around the Bow valley, west of Calgary, Canada. This area is approximate 90,000 km^2 and the mesh has approximately 130k triangles. The mesh is shown with triangle edges coloured in grey, and the faces coloured by elevation. Areas of significant heterogeneity, such as the steep ridge lines and the tree-line, have increased triangle density. Areas with low heterogeneity, such as the valleys, have large triangles. This mesh is currently used for the [SnowCast](http://www.snowcast.ca) product.

![](images/mesh.png)


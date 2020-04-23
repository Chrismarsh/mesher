Output
======

Mesher creates a directory with the same name as the configuration file name or a user-defined name as configured by ``user_output_dir`` (see :ref:`configuration:Outputs`). Within this folder is a subdirectory with the same name as that of the input DEM.

In the primarily output directory are three JSON format files:
- ``*.ic`` that contains the initial conditions
- ``*.mesh`` that contains the primary mesh
- ``*.param`` that contains the per-triangle parameter values

If no parameter or initial conditions are given, then these files are empty.


In the second directory are the reprojected files (```*_projected```), intermediary files, and the triangulation shape file (```*_USM.shp```). A ```*.vtu``` file is also created for visualizing in 3D in `Paraview <http://www.paraview.org>`_.

.. note::
   All indicies are 0-indexed!

.mesh
******

The primary mesh file's structure is as follows:
::

 {
   "mesh": { 
            "vertex": [ [double,double,double], [double,double,double], ... ],
            "nvertex": int,
            "neigh": [ [int,int,int], [int,int,int], ... ],
            "elem": [ [int,int,int], [int,int,int], ... ],
            },
   "is_geographic": int,
   "proj4": string,
   "UTM_zone": int,
   "nelem": int
 }

.. confval:: mesh.vertex

   :type: list of [double,double,double] 

List of triangle vertexes. Triangles may share a vertex, thus it isn't stored twice. A vertex has [x,y,z] coordinates in the output coordinate system.


.. confval:: mesh.nvertex

   :type: int

Number of verticies.

.. confval:: mesh.neigh

   :type: list of [int,int,int]

The neighbours to each triangle. Each int is an index into the ``elem`` list. 

A value of ``-1`` means "missing". 

.. confval:: mesh.elem

   :type: list of [int,int,int]

Verticies are connected via edges to form triangle elements. These indicies are indexes into the vertex array. 


For example, consider the :ref:`examples:flat` mesh. It has 2 triangles in the domain. The ``.mesh`` file looks like:
::

   {
       "mesh": {
           "vertex": [
               [
                   488479.5,
                   6713528.0,
                   1000.0
               ],
               [
                   490527.5,
                   6713528.0,
                   1000.0
               ],
               [
                   490527.5,
                   6711480.0,
                   1000.0
               ],
               [
                   488479.5,
                   6711480.0,
                   1000.0
               ]
           ],
           "nvertex": 4,
           "neigh": [
               [
                   1,
                   -1,
                   -1
               ],
               [
                   -1,
                   0,
                   -1
               ]
           ],
           "elem": [
               [
                   1,
                   0,
                   2
               ],
               [
                   0,
                   3,
                   2
               ]
           ],
           "is_geographic": 0,
           "proj4": "+proj=utm +zone=8 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ",
           "UTM_zone": 8,
           "nelem": 2
       }
   }

For example, triangle0 is give by mesh.elem[0]:
::

   [
       1,
       0,
       2
   ]

This says that to define triangle0, use the 1st, the 0th, and 2nd vertex elements. That is:
::

   triangle0 = [490527.5, 6713528.0, 1000.0], [488479.5, 6713528.0, 1000.0], [490527.5, 6711480.0, 1000.0]

It has a neighbour ``mesh.neigh[0]`` 
::

   [
      1,
      -1,
      -1
   ]

which says that edge0 is the only edge with a neighbour, and it is triangle1.



.param
*******

The parameter file's schema is:
::

   {
       "param1": [
           double/int,
           double/int,,
           ...
       ],
      "param2": [
           double/int,
           double/int,,
           ...
       ]
   }


.. confval:: param1

   :type: List of length ``length(elem)``

Each item of the list at index *i* corresponds to a value for that parameter for triangle_*i*.


For example, that parameter file for ideal is
::

   {
       "area": [
           2097152.0,
           2097152.0
       ]
   }


Thus the area of triangle0 is 2097152 m^2 or (1448m)^2. The latter is more easily used for approximate comparisons against raster cell sizes.




.ic
****

The initial conditions file's schema is:
::

   {
       "ic1": [
           double/int,
           double/int,,
           ...
       ],
      "ic2": [
           double/int,
           double/int,,
           ...
       ]
   }


.. confval:: ic1

   :type: List of length ``length(elem)``

Each item of the list at index *i* corresponds to a value for that initial condition for triangle_*i*.


.vtu
****
Within the intermediary folder is a ``.vtu`` file. This can be opened in Paraview to provide feedback on how the mesh generation went, and if tweaks should be done.


.shp
*****
Within the intermediary folder is an ESRI ``.shp`` file with the ``_USM.shp`` suffix. This can be loaded directly into a GIS tool.





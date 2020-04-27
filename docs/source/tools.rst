Tools
=======

mesh2vtu.py
***********

Combines and converts a ``.mesh`` file and a ``.param`` file into a corresponding ``.vtu`` file

Usage:
::

   mesh2vtu.py input.mesh input.param output.vtu


meshmerge.py
*************

Combines multiple parameter files (``.param``) into one. All the files are validated to ensure each parameter has the same number of elements.

Usage:
::

   meshmerge.py [-h] [-o outfile] infile [infile ...]

.. confval:: infile

   :type: string

Path to the *n* input parameter files to merge into one

.. confval:: outfile

   :type: string

Optional path to the output parameter file. If not given, defaults to a filename of concatenated input files.


mesherpermuation.py
*******************
Because of how the triangles are generated and subsequently written to file, triangles that are close in space may not be close in the file. Therefore, if the mesh is used in a numerical model, it may result in inefficient access. This tool computes the permutation of a mesh that minimizes the bandwidth of the connectivity matrix. 


This tool adds a new field ``cell_global_id`` to the ``.mesh`` file. A consuming program can sort the triangles (i.e., ``mesh.elem`` see: :ref:`output:.mesh`) by these ids. Then, spatially close triangles will also be close in memory. This may be visualized by plotting the nearest neighbour connectivity before and after reodering.

**Before optimization**

.. image:: images/orig_order.png

**After optimization**

.. image:: images/post_order.png





Usage:
::

   mesherpermuation.py [-h] [-o outfile] [-t type] [-i infile]

.. confval:: outfile

File for output. Overwrite input file


.. confval:: type


``rcm`` Performing RCM bandwidth minimization (default)

``nd`` Performing ND fill-in minimization
























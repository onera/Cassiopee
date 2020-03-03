.. OCC documentation master file


OCC: interface to open-cascade
===============================

Preamble
########

OCC performs reading of IGES or STEP files with open-cascade. It returns
a triangular mesh.

To use the module with the Converter.array interface::

   import OCC

To use the module with the CGNS/Python interface::

    import OCC.PyTree as OCC


.. py:module:: OCC


List of functions
##################

**-- CAD/mesh conversion**

.. autosummary::

    OCC.convertCAD2Arrays
    OCC.PyTree.convertCAD2PyTree

   

Contents
#########


CAD/mesh conversion
----------------------------


.. py:function:: OCC.convertCAD2Arrays(fileName, format='fmt_iges', h=0., chordal_err=0., growth_ratio=0., algo=0)

    Read a CAD and return arrays.

       
    :param fileName: CAD file name
    :type fileName: string
    :param format: file format ('fmt_iges' or 'fmt_step')
    :type format: string
    :param h: step size on output mesh. If 0., automatic setting [algo=0].
    :type h: float
    :param chordal_error: max error between CAD and mesh. Result in curvature adaptation. If 0., automatic setting.
    :type chordal_error: float
    :param growth_ratio: max growth ratio between adjacent triangles [algo=0].
    :type growth_ratio: float
    :param algo: algo=0: mesh with regular triangles, algo=1: mesh with only respect to curvature.
    :type algo: int
    :rtype: a list of TRI arrays

    *Example of use:*

    * `Read a CAD (array) <Examples/OCC/convertCAD2Arrays.py>`_:

    .. literalinclude:: ../build/Examples/OCC/convertCAD2Arrays.py

---------------------------------------

.. py:function:: OCC.PyTree.convertCAD2PyTree(fileName, format='fmt_iges', h=0., chordal_err=0., growth_ratio=0., algo=0)

    Read a CAD and return zone.

    :param fileName: CAD file name
    :type fileName: string
    :param format: file format ('fmt_iges' or 'fmt_step')
    :type format: string
    :param h: step size on output mesh. If 0., automatic setting [algo=0].
    :type h: float
    :param chordal_error: max error between CAD and mesh. Result in curvature adaptation. If 0., automatic setting.
    :type chordal_error: float
    :param growth_ratio: max growth ratio between adjacent triangles [algo=0].
    :type growth_ratio: float
    :param algo: algo=0: mesh with regular triangles, algo=1: mesh with only respect to curvature.
    :type algo: int
    :rtype: CGNS pyTree

    *Example of use:*

    * `Read a CAD (pyTree) <Examples/OCC/convertCAD2PyTreePT.py>`_:

    .. literalinclude:: ../build/Examples/OCC/convertCAD2PyTreePT.py

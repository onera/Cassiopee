.. OCC documentation master file

:tocdepth: 2


OCC: interface to open-cascade
===============================

Preamble
########

OCC performs CAD operations and surface meshing using open-cascade. 

To use the module with the Converter.array interface::

   import OCC

To use the module with the CGNS/Python interface::

    import OCC.PyTree as OCC


.. py:module:: OCC


List of functions
##################

**-- CAD to surface mesh**

.. autosummary::
   :nosignatures:

    OCC.PyTree/meshAll

**-- CAD to surface mesh (legacy functions)**

.. autosummary::
   :nosignatures:

    OCC.convertCAD2Arrays
    OCC.PyTree.convertCAD2PyTree

**-- CAD manipulation**

.. autosummary::
   :nosignatures:

    OCC.readCAD
    OCC.writeCAD
    OCC.getNbEdges
    OCC.getNbFaces
    OCC.getFaceArea
    OCC._translate
    OCC._rotate
    OCC._splitFaces
    OCC._mergeFaces


Contents
#########


CAD to surface mesh (legacy functions)
---------------------------------------

.. py:function:: OCC.convertCAD2Arrays(fileName, format='fmt_iges', h=0., chordal_err=0., growth_ratio=0., algo=1)

    Read a CAD and return arrays.

    :param fileName: CAD file name
    :type fileName: string
    :param format: file format ('fmt_iges' or 'fmt_step')
    :type format: string
    :param h: step size on output mesh. If 0., automatic setting [algo=1,2].
    :type h: float
    :param chordal_error: max error between CAD and mesh. Result in curvature adaptation. If 0., automatic setting.
    :type chordal_error: float
    :param growth_ratio: max growth ratio between adjacent triangles [algo=1,2].
    :type growth_ratio: float
    :param algo: algo=0: mesh with only respect to curvature, algo=1 or algo=2: mesh with regular triangles.
    :type algo: int
    :rtype: a list of TRI arrays

    *Example of use:*

    * `Read a CAD (array) <Examples/OCC/convertCAD2Arrays.py>`_:

    .. literalinclude:: ../build/Examples/OCC/convertCAD2Arrays.py

---------------------------------------

.. py:function:: OCC.PyTree.convertCAD2PyTree(fileName, format='fmt_iges', h=0., chordal_err=0., growth_ratio=0., algo=1)

    Read a CAD and return a zone.

    :param fileName: CAD file name
    :type fileName: string
    :param format: file format ('fmt_iges' or 'fmt_step')
    :type format: string
    :param h: step size on output mesh. If 0., automatic setting [algo=1,2].
    :type h: float
    :param chordal_error: max error between CAD and mesh. Result in curvature adaptation. If 0., automatic setting.
    :type chordal_error: float
    :param growth_ratio: max growth ratio between adjacent triangles [algo=1,2].
    :type growth_ratio: float
    :param algo: algo=0: mesh with only respect to curvature, algo=1 or algo=2: mesh with regular triangles. 
    :type algo: int
    :rtype: CGNS pyTree

    *Example of use:*

    * `Read a CAD (pyTree) <Examples/OCC/convertCAD2PyTreePT.py>`_:

    .. literalinclude:: ../build/Examples/OCC/convertCAD2PyTreePT.py


CAD manipulation
----------------------

.. py:function:: OCC.readCAD(fileName, format='fmt_step')

    Read a CAD file and return a CAD hook.

    :param fileName: CAD file name
    :type fileName: string
    :param format: file format ('fmt_iges' or 'fmt_step')
    :type format: string
    :rtype: CAD hook

    *Example of use:*

    * `Read a CAD <Examples/OCC/readCAD.py>`_:

    .. literalinclude:: ../build/Examples/OCC/readCAD.py

------------------------------------------

.. py:function:: OCC.writeCAD(hook, fileName, format='fmt_step')

    Write a CAD hook to a file.

    :param hook: CAD hook
    :type hook: CAD hook
    :param fileName: CAD file name
    :type fileName: string
    :param format: file format ('fmt_iges' or 'fmt_step')
    :type format: string

    *Example of use:*

    * `Write a CAD <Examples/OCC/writeCAD.py>`_:

    .. literalinclude:: ../build/Examples/OCC/writeCAD.py

------------------------------------------

.. py:function:: OCC.getNbEdges(hook)

    Return the number of edges in a CAD hook.

    :param hook: CAD hook
    :type hook: CAD hook
    :rtype: int

    *Example of use:*

    * `Get the number of edges <Examples/OCC/getNbEdges.py>`_:

    .. literalinclude:: ../build/Examples/OCC/getNbEdges.py

------------------------------------------

.. py:function:: OCC.getNbFaces(hook)

    Return the number of faces in a CAD hook.

    :param hook: CAD hook
    :type hook: CAD hook
    :rtype: int

    *Example of use:*

    * `Get the number of faces <Examples/OCC/getNbFaces.py>`_:

    .. literalinclude:: ../build/Examples/OCC/getNbFaces.py

------------------------------------------

.. py:function:: OCC._translate(hook, vector, listOfFaces=None)

    Translate a CAD hook by a given vector.

    :param hook: CAD hook
    :type hook: CAD hook
    :param vector: translation vector (dx, dy, dz)
    :type vector: tuple of floats
    :param listOfFaces: if None, translate all else translate only given faces
    :type listOfFaces: list of face indices (starts 1)

    *Example of use:*

    * `Translate a CAD <Examples/OCC/translate.py>`_:

    .. literalinclude:: ../build/Examples/OCC/translate.py

------------------------------------------

.. py:function:: OCC._rotate(hook, Xc, axis, angle, listOfFaces=None)

    Rotate a CAD hook around a given axis by a given angle.

    :param hook: CAD hook
    :type hook: CAD hook
    :param Xc: rotation center (x, y, z)
    :type Xc: tuple of floats
    :param axis: rotation axis
    :type axis: tuple of floats
    :param angle: rotation angle in degrees
    :type angle: float
    :param listOfFaces: if None, rotate all else rotate only given faces
    :type listOfFaces: list of face indices (starts 1)

    *Example of use:*

    * `Rotate a CAD <Examples/OCC/rotate.py>`_:

    .. literalinclude:: ../build/Examples/OCC/rotate.py

------------------------------------------

.. py:function:: OCC._splitFaces(hook, area)

    Split faces in a CAD hook.

    :param hook: CAD hook
    :type hook: CAD hook
    :param area: split each face if area greater than this value
    :type area: float

    *Example of use:*

    * `Split faces <Examples/OCC/splitFaces.py>`_:

    .. literalinclude:: ../build/Examples/OCC/splitFaces.py

------------------------------------------

.. py:function:: OCC._mergeFaces(hook, listFaces=None)

    Merge faces in a CAD hook.

    :param hook: CAD hook
    :type hook: CAD hook
    :param listFaces: if None, merge all faces else list of faces indices to merge
    :type listFaces: list of face indices (starts 1)

    *Example of use:*

    * `Merge faces <Examples/OCC/mergeFaces.py>`_:

    .. literalinclude:: ../build/Examples/OCC/mergeFaces.py

------------------------------------------

.. py:function:: OCC.getFaceArea(hook, listFaces=[])

    Return the area of given faces.

    :param hook: CAD hook
    :type hook: CAD hook
    :param listFaces: list of faces number to calculate the area
    :type listFaces: list of face index (starts 1)
    :rtype: float

    *Example of use:*

    * `Get face area <Examples/OCC/getFaceArea.py>`_:

    .. literalinclude:: ../build/Examples/OCC/getFaceArea.py


    



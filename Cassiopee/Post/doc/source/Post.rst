.. Post documentation master file

:tocdepth: 2


Post: post-processing module
=========================================

Preamble
########

This module provides post-processing tools for  CFD simulations.
It manipulates arrays (as defined in Converter documentation)
or CGNS/Python trees (or pyTree, as defined in Converter/Internal documentation)
as data structures.

This module is part of Cassiopee, a free open-source
pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Post module::

   import Post

For use with the pyTree interface::

    import Post.PyTree as Post


.. py:module:: Post

List of functions
##################

**-- Modifying/creating Variables**

.. autosummary::
   :nosignatures:

    Post.renameVars
    Post.PyTree.importVariables
    Post.computeVariables
    Post.computeExtraVariable
    Post.PyTree.computeWallShearStress
    Post.computeGrad
    Post.computeGrad2
    Post.computeGradLSQ
    Post.computeNormGrad
    Post.computeDiv
    Post.computeDiv2
    Post.computeCurl
    Post.computeNormCurl
    Post.computeDiff


**-- Solution selection**

.. autosummary::
   :nosignatures:

    Post.selectCells
    Post.selectCells2
    Post.exteriorVertices
    Post.interiorFaces
    Post.exteriorFaces
    Post.exteriorFacesStructured
    Post.exteriorElts
    Post.frontFaces
    Post.sharpEdges
    Post.silhouette
    Post.coarsen
    Post.refine
    Post.computeIndicatorValue
    Post.computeIndicatorField


**-- Solution extraction**

.. autosummary::
   :nosignatures:

    Post.extractPoint
    Post.extractPlane
    Post.extractMesh
    Post.projectCloudSolution
    Post.zipper
    Post.usurp

**-- Streams/Isos**

.. autosummary::
   :nosignatures:

    Post.streamLine
    Post.streamRibbon
    Post.streamSurf
    Post.isoLine
    Post.isoSurf
    Post.isoSurfMC

**-- Solution integration**

.. autosummary::
   :nosignatures:

    Post.integ
    Post.integNorm
    Post.integNormProduct
    Post.integMoment
    Post.integMomentNorm


Contents
#########

Modifying/creating variables
------------------------------


.. py:function:: Post.renameVars(t, oldVarNameList, newVarNameList)

    Rename a list of variables with new variable names.
    Exists also as in place function (_renameVars) that modifies t and returns None.

    :param t:  Input data
    :type  t:  [array, arrays] or [zone, list of zones, base, tree]
    :param oldVarNameList: list of variables to rename
    :type  oldVarNameList: list of strings
    :param newVarNameList: list of new variable names
    :type  newVarNameList: list of strings
    :return: reference copy of input
    :rtype: identical to t

    *Example of use:*

    * `Rename variables (array) <Examples/Post/renameVars.py>`_:

    .. literalinclude:: ../build/Examples/Post/renameVars.py

    * `Rename variables (pyTree) <Examples/Post/renameVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/renameVarsPT.py


---------------------------------------

.. py:function:: Post.PyTree.importVariables(t1, t2, method=0, eps=1.e-6, addExtra=1)

    Variables located at nodes and/or centers can be imported from a pyTree t1
    to a pyTree t2.
    If one variable already exists in t2, it is replaced by the same
    variable from t1.
    If method=0, zone are matched from names, if method=1, zones are
    matched from coordinates with a tolerance eps, if method=2, zones
    are taken in the given order of t1 and t2 (must match one by one).
    If addExtra=1, unmatched zones are added to a base named 'EXTRA'.

    :param t1:  Input data
    :type  t1:  pyTree
    :param t2:  Input data
    :type  t2:  pyTree
    :return: reference copy of t2
    :rtype: pyTree

    *Example of use:*

    * `Import variables to a tree (pyTree) <Examples/Post/importVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/importVariablesPT.py

---------------------------------------

.. py:function:: Post.computeVariables(a, varList, gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6, Cs=110.4, mus=1.76e-5, Ts=273.15)

    New variables can be computed from conservative variables.
    The list of the names of the variables to compute must be provided.
    The computation of some variables (e.g. viscosity) require some constants as input data.
    In the pyTree version, if a reference state node is defined in the pyTree, then the corresponding reference
    constants are used. Otherwise, they must be specified as an argument of the function.
    Exists also as in place version (_computeVariables) that modifies a and returns None.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varList: list of variable names (can be preceded by 'nodes:' or 'centers:')
    :type varList: list of strings
    :rtype:  identical to input

    The constants are:

    - 'gamma' for the specific heat ratio: :math:`\gamma`;
    - 'rgp' for the perfect gas constant: :math:`R = (\gamma-1) \times C_v`;
    - 'betas' and 'Cs' (Sutherland's law constants), or 'Cs','Ts' and 'mus';
    - 's0' for a constant entropy, defined by: :math:`s_0 = s_{ref} - R \frac{\gamma}{\gamma-1} ln(T_{ref}) + R\ ln(P_{ref})` where :math:`\ s_{ref}, T_{ref}` and :math:`P_{ref}` are defined for a reference state.

    Computed variables are defined by their CGNS names:

    * 'VelocityX', 'VelocityY', 'VelocityZ' for components of the absolute velocity,
    * 'VelocityMagnitude' for the absolute velocity magnitude,
    * 'Pressure' for the static pressure (requires: gamma),
    * 'Temperature' for the static temperature (requires: gamma, rgp),
    * 'Enthalpy' for the enthalpy (requires: gamma),
    * 'Entropy' for the entropy (requires: gamma, rgp, s0),
    * 'Mach' for the Mach number (requires: gamma),
    * 'ViscosityMolecular' for the fluid molecular viscosity (requires: gamma, rgp, Ts, mus, Cs),
    * 'PressureStagnation' for stagnation pressure(requires: gamma),
    * 'TemperatureStagnation' for stagnation temperature (requires: gamma, rgp),
    * 'PressureDynamic' for dynamic pressure (requires: gamma).


    *Example of use:*

    * `Compute variables (array) <Examples/Post/computeVariables.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeVariables.py

    .. note:: In the pyTree version, if the variable name is prefixed by 'centers:' then the variable is computed at centers only (e.g. 'centers:Pressure'), and if it is not prefixed, then the variable is computed at nodes.

    * `Compute variables (pyTree) <Examples/Post/computeVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeVariablesPT.py

---------------------------------------

.. py:function:: Post.computeExtraVariable(a, varName, gamma=1.4, rgp=287.053, Cs=110.4, mus=1.76e-5, Ts=273.15)

    Compute more advanced variables from conservative variables.
    'varName' can be:

    - Vorticity,
    - VorticityMagnitude,
    - QCriterion,
    - ShearStress,
    - SkinFriction,
    - SkinFrictionTangential

    The computation of the shear stress requires  gamma, rgp, Ts, mus, Cs as input data.
    In the pyTree version, if a reference state node is defined in the pyTree, then thecorresponding reference
    constants are used. Otherwise, they must be specified as an argument of the function.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varName: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varName: string
    :rtype: identical to input

    *Example of use:*

    * `Extra variables computation (array) <Examples/Post/computeExtraVariable.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeExtraVariable.py

    * `Extra variables computation (pyTree) <Examples/Post/computeExtraVariablePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeExtraVariablePT.py

---------------------------------------

.. py:function:: Post.PyTree.computeWallShearStress(t)

    Compute the shear stress at wall boundaries provided the velocity gradient is already computed.
    The problem dimension and the reference state must be provided in t, defining the skin mesh.

    Exists also as in place version (_computeWallShearStress) that modifies t and returns None.

    The function is only available in the pyTree version.

    :param t: input data
    :type  t: pyTree, base, zone, list of zones
    :rtype:  identical to input

    *Example of use:*

    * `Wall shear stress computation (pyTree) <Examples/Post/computeWallShearStressPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeWallShearStressPT.py

---------------------------------------

.. py:function:: Post.computeGrad(a, varname)

    Compute the gradient (:math:`\nabla x, \nabla y, \nabla z`) of a field of name *varname*
    defined in *a*. The returned field is located at cell centers.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varname: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varname: string
    :rtype:  identical to input

    *Example of use:*

    * `Gradient of density field (array) <Examples/Post/computeGrad.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeGrad.py

    * `Gradient of density field (pyTree) <Examples/Post/computeGradPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeGradPT.py

---------------------------------------

.. py:function:: Post.computeGrad2(a, varname)

    Compute the gradient (:math:`\nabla x, \nabla y, \nabla z`) at cell centers for a field of name *varname* located at cell centers.

    Using Converter.array interface:
    ::

        P.computeGrad2(a, ac, indices=None, BCField=None)

    *a* denotes the mesh, *ac* denotes the fields located at centers.
    indices is a numpy 1D-array of face list, BCField is the corresponding numpy array of face fields. They are used to force a value at some faces before computing the gradient.

    Using the pyTree version:
    ::

        P.computeGrad2(a, varname)

    The variable name must be located at cell centers.
    Indices and BCFields are automatically extracted from BCDataSet nodes:
    if a BCDataSet node is defined for a BC of the pyTree, the corresponding face fields
    are imposed when computing the gradient.
    If volume has already been computed and volume field is present in tree, it is not recomputed for the gradient computation (only NGON cases up to now).

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varname: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varname: string
    :rtype:  identical to input

    *Example of use:*

    * `Gradient of density field with computeGrad2 (array) <Examples/Post/computeGrad2.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeGrad2.py

    * `Gradient of density field with computeGrad2 (pyTree) <Examples/Post/computeGradPT2.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeGrad2PT.py

---------------------------------------

.. py:function:: Post.computeGradLSQ(a, varNames)

    Compute the gradient (:math:`\nabla x, \nabla y, \nabla z`) at cell centers for a list of fields of names *varNames* located at cell centers. Only supports NGon meshes.

    Using the pyTree version:
    ::

        a = P.computeGradLSQ(a, varNames)

    :param a: Input data
    :type  a: [pyTree, base, zone, list of zones]
    :param varname: list of variable names
    :type varname: string
    :rtype: identical to input

    *Example of use:*

    * `Gradient of fields f and g defined by functions F and G with computeGradLSQ (pyTree) <Examples/Post/computeGradLSQPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeGradLSQPT.py

---------------------------------------

.. py:function:: Post.computeDiv(a, varname)

    Compute the divergence :math:`\nabla\cdot\left(\vec{\bullet}\right)` of a field defined by its
    component names ['vectX','vectY','vectZ'] defined in *a*.
    The returned field is located at cell centers.

    Using Converter.array interface:
    ::

        P.computeDiv(a, ['vectX','vectY','vectZ'])

    Using the pyTree version:
    ::

        P.computeDiv(a, 'vect')

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varname: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varname: string
    :rtype:  identical to input

    *Example of use:*

    * `Divergence of a vector field (array)  with computeDiv <Examples/Post/computeDiv.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeDiv.py

    * `Divergence of a vector field (pyTree)  with computeDiv <Examples/Post/computeDivPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeDivPT.py

---------------------------------------

.. py:function:: Post.computeDiv2(a, varname)

    compute the divergence :math:`\nabla\cdot\left(\vec{\bullet}\right)` at cell
    centers for a vector field defined by its variable names ['vectX','vectY','vectZ']
    located at cell centers.

    Using Converter.array interface:
    ::

        P.computeDiv2(a, ac, indices=None, BCField=None)

    *a* denotes the mesh, *ac* denotes the components of the vector field located at centers.
    indices is a numpy 1D-array of face list, BCField is the corresponding numpy array of face fields.
    They are used to force a value at some faces before computing the gradients.

    Using the pyTree version:
    ::

        P.computeDiv2(a, 'vect')
        P.computeDiv2(a, ['vect1', 'vect2'])

    The variable name must be located at cell centers.
    Indices and BCFields are automatically extracted from BCDataSet nodes:
    if a BCDataSet node is defined for a BC of the pyTree, the corresponding face fields
    are imposed when computing the gradient.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varname: variable name(s) (can be preceded by 'nodes:' or 'centers:')
    :type varname: [string, list of strings]
    :rtype:  identical to input

    *Example of use:*

    * `Divergence of a vector field (array) with computeDiv2 <Examples/Post/computeDiv2.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeDiv2.py

    * `Divergence of a vector field (pyTree) with computeDiv2 <Examples/Post/computeDiv2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeDiv2PT.py

---------------------------------------

.. py:function:: Post.computeNormGrad(a, varname)

    Compute the norm of gradient (:math:`\nabla x, \nabla y, \nabla z`) of a field of name varname defined in a. The returned field 'grad'+varname and is located at cell centers. **(???)**

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varname: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varname: string
    :rtype:  identical to input

    *Example of use:*

    * `Norm of gradient of density (array) <Examples/Post/computeNormGrad.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeNormGrad.py

    * `Norm of gradient of density (pyTree) <Examples/Post/computeNormGradPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeNormGradPT.py

---------------------------------------

.. py:function:: Post.computeCurl(a, ['vectx','vecty','vectz'])

    Compute curl of a 3D vector defined by its variable names
    ['vectx','vecty','vectz'] in a.
    The returned field is defined at cell centers for structured grids and elements centers for unstructured grids.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param vect*: variable name defining the 3D vector
    :type vect*: string
    :rtype:  identical to input


    *Example of use:*

    * `Curl of momentum field (array) <Examples/Post/computeCurl.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeCurl.py

    * `Curl of momentum field (pyTree) <Examples/Post/computeCurlPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeCurlPT.py

---------------------------------------

.. py:function:: Post.computeNormCurl(a, ['vectx','vecty','vectz'])

    Compute the norm of the curl of a 3D vector defined by its variable names
    ['vectx','vecty','vectz'] in a.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param vect*: variable name defining the 3D vector
    :type vect*: string
    :rtype:  identical to input

    *Example of use:*

    * `Norm of the curl of momentum field (array) <Examples/Post/computeNormCurl.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeNormCurl.py

    * `Norm of the curl of momentum field (pyTree) <Examples/Post/computeNormCurlPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeNormCurlPT.py

---------------------------------------

.. py:function:: Post.computeDiff(a, varname)

    Compute the difference between neighbouring cells of a scalar field defined by its variable varname in a.
    The maximum of the absolute difference among all directions is kept.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varname: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varname: string
    :rtype:  identical to input

    *Example of use:*

    * `Difference of density field (array) <Examples/Post/computeDiff.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeDiff.py

    * `Difference of density field  (pyTree) <Examples/Post/computeDiffPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeDiffPT.py

---------------------------------------

Solution selection
-------------------

.. py:function:: Post.selectCells(a, F, ['var1', 'var2'], strict=0, cleanConnectivity=True)

    Select cells with respect to a given criterion.
    If strict=0, the cell is selected if at least one of the cell vertices satisfies the criterion.
    If strict=1, the cell is selected if all the cell vertices satisfy the criterion.
    The criterion can be defined as a python function returning True (=selected) or False (=not selected):
    ::

        P.selectCells(a, F, ['var1', 'var2'], strict=0)

    or by a formula:
    ::

        P.selectCells(a, '{x}+{y}>2', strict=0)

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param F: cells selection criterion
    :type F: function
    :param var*: arguments of function F
    :type var*: string
    :param strict: selection mode (0 or 1)
    :type strict: integer
    :param cleanConnectivity: if True, connectivity is cleaned
    :type cleanConnectivity: True or False
    :rtype: identical to input

    *Example of use:*

    * `Cell selection in a mesh (array) <Examples/Post/selectCells.py>`_:

    .. literalinclude:: ../build/Examples/Post/selectCells.py

    * `Cell selection in a mesh  (pyTree) <Examples/Post/selectCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/selectCellsPT.py

---------------------------------------

.. py:function:: Post.selectCells2(a, tag, strict=0)

    Select cells according to a field defined by a variable 'tag' (=1 if selected, =0 if not selected).
    If 'tag' is located at centers, only cells of tag=1 are selected.
    If 'tag' is located at nodes and 'strict'=0, the cell is selected if at least one of the cell vertices is tag=1.
    If 'tag' is located at nodes and 'strict'=1, the cell is selected if all the cell vertices is tag=1.
    In the array version, the tag is an array. In the pyTree version, the tag must be defined in a 'FlowSolution_t' type node
    located at cell centers or nodes.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tag: variable name
    :type tag: string
    :param strict: selection mode (0 or 1)
    :type strict: integer
    :param cleanConnectivity: if True, connectivity is cleaned
    :type cleanConnectivity: True or False

    :rtype: identical to input

    *Example of use:*

    * `Cell selection in a mesh with selectCells2 (array) <Examples/Post/selectCells2.py>`_:

    .. literalinclude:: ../build/Examples/Post/selectCells2.py

    * `Cell selection in a mesh with selectCells 2 (pyTree) <Examples/Post/selectCells2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/selectCells2PT.py

-----------------------------------------

.. py:function:: Post.exteriorVertices(a, indices=None)

    Select the exterior vertices of a mesh, and return them in a single 
    unstructured NODE zone. If indices=[], the indices of the original exterior 
    vertices are returned.

    Exists also as in place version (_exteriorVertices) that modifies a and return None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param indices: indices of original exterior vertices
    :type indices: list of integers
    :rtype: zone

    *Example of use:*

    * `Select exterior vertices (array) <Examples/Post/exteriorVertices.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorVertices.py

    * `Select exterior vertices (pyTree) <Examples/Post/exteriorVerticesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorVerticesPT.py
        
---------------------------------------

.. py:function:: Post.interiorFaces(a, strict=0)

    Select the interior faces of a mesh. Interior faces are faces with
    two neighbouring elements. If 'strict' is set to 1, select the interior faces
    that have only interior nodes.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param strict: selection mode (0 or 1)
    :type strict: integer
    :rtype: identical to input

    *Example of use:*

    * `Select interior faces (array) <Examples/Post/interiorFaces.py>`_:

    .. literalinclude:: ../build/Examples/Post/interiorFaces.py

    * `Select interior faces (pyTree) <Examples/Post/interiorFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/interiorFacesPT.py

-----------------------------------------

.. py:function:: Post.exteriorFaces(a, indices=None)

    Select the exterior faces of a mesh, and return them in a single unstructured zone. If indices=[], the
    indices of the original exterior faces are returned.
    For structured grids, indices are the global index containing i faces, then j faces, then k faces, starting from 0.
    For NGON grids, indices are the NGON face indices, starting from 1.

    Exists also as in place version (_exteriorFaces) that modifies a and return None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param indices: indices of original exterior faces
    :type indices: list of integers
    :rtype: zone

    *Example of use:*

    * `Select exterior faces (array) <Examples/Post/exteriorFaces.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorFaces.py

    * `Select exterior faces (pyTree) <Examples/Post/exteriorFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorFacesPT.py

---------------------------------------

.. py:function:: Post.exteriorFacesStructured(a)

    Select the exterior faces of a structured mesh as a list of structured meshes.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: zone

    *Example of use:*

    * `Select structured exterior faces (array) <Examples/Post/exteriorFacesStructured.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorFacesStructured.py

    * `Select structured exterior faces (pyTree) <Examples/Post/exteriorFacesStructuredPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorFacesStructuredPT.py

---------------------------------------

.. py:function:: Post.exteriorElts(a)

    Select the exterior elements of a mesh, that is the first border fringe of cells.

    Exists also as in place version (_exteriorElts) that modifies a and return None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: identical to input

    *Example of use:*

    * `Select exterior elements (array) <Examples/Post/exteriorElts.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorElts.py

    * `Select exterior elements (pyTree) <Examples/Post/exteriorEltsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exteriorEltsPT.py

---------------------------------------

.. py:function:: Post.frontFaces(a, tag)

    Select faces that are located at the boundary where a tag indicator change from 0 to 1.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tag: variable name
    :type tag: string
    :rtype: zone

    *Example of use:*

    * `Select a front in a tag (array) <Examples/Post/frontFaces.py>`_:

    .. literalinclude:: ../build/Examples/Post/frontFaces.py

    * `Select a front in a tag (pyTree) <Examples/Post/frontFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/frontFacesPT.py

---------------------------------------

.. py:function:: Post.sharpEdges(A, alphaRef=30.)

    Return sharp edges arrays starting from surfaces or contours.
    Adjacent cells having an angle deviating from more than alphaRef to 180 degrees are considered as sharp.

    :param A: input data
    :type A: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param alphaRef: split angle
    :type alphaRef: float
    :rtype: list of arrays / zones **??**

    *Example of use:*

    * `Detect sharp edges of a surface (array) <Examples/Post/sharpEdges.py>`_:

    .. literalinclude:: ../build/Examples/Post/sharpEdges.py

    * `Detect sharp edges of a surface (pyTree) <Examples/Post/sharpEdgesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/sharpEdgesPT.py

---------------------------------------

.. py:function:: Post.silhouette(A, vector=[1.,0.,0.])

    Return silhouette arrays starting from surfaces or contours, according to a direction vector.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param vector: direction vector
    :type vector: 3-tuple of floats
    :rtype: identical to input

    *Example of use:*

    * `Detect silhouette of a surface (array) <Examples/Post/silhouette.py>`_:

    .. literalinclude:: ../build/Examples/Post/silhouette.py

    * `Detect silhouette of a surface (pyTree) <Examples/Post/silhouettePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/silhouettePT.py

---------------------------------------

.. py:function:: Post.coarsen(a, indicName='indic', argqual=0.25, tol=1.e-6)

    Coarsen a triangle mesh by providing a coarsening indicator, which is 1 if the element must be coarsened, 0 elsewhere.
    Triangles are merged by edge contraction, if tagged to be coarsened
    by indic and if new triangles deviate less than tol to the original triangle.
    Required mesh quality is controled by argqual: argqual equal to 0.5
    corresponds to an equilateral triangle,
    whereas a value near zero corresponds to a bad triangle shape.

    **Array version**: an indic i-array must be provided, whose dimension ni is equal to the number of elements in the initial triangulation:
    ::

        b = P.coarsen(a, indic, argqual=0.1, tol=1.e6)


    :param a: input data
    :type a: array, list of arrays
    :param indic: tagged element (0 or 1)
    :type indic: i-array
    :rtype: identical to input

    **PyTree version**: indic is stored as a solution located at centers:
    ::

        b = P.coarsen(a, indicName='indic', argqual=0.25, tol=1.e-6)

    :param a: input data
    :type a: pyTree, base, zone, list of zones
    :param indicName: tag variable name
    :type indicName: string
    :rtype: identical to input

    *Example of use:*

    * `Coarsen all cells in a 2D mesh (array) <Examples/Post/coarsen.py>`_:

    .. literalinclude:: ../build/Examples/Post/coarsen.py

    * `Coarsen all cells in a 2D mesh (pyTree) <Examples/Post/coarsenPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/coarsenPT.py

---------------------------------------

.. py:function:: Post.refine

    Refine a triangle mesh by providing a refinement indicator, which is 1 if the element must be refined, 0 elsewhere.

    **Array version**: an indic i-array must be provided, whose dimension ni
    is equal to the number of elements in the initial triangulation:
    ::

        b = P.refine(a, indic)

    **PyTree version**: indic is stored as a solution located at centers:
    ::

        b = P.refine(a, indicName='indic')

    *Example of use:*

    * `Refine all cells in a 2D mesh (array) <Examples/Post/refine.py>`_:

    .. literalinclude:: ../build/Examples/Post/refine.py

    * `Refine all cells in a 2D mesh (pyTree) <Examples/Post/refinePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/refinePT.py

---------------------------------------

.. py:function:: Post.refine(a, w=1./64.)

    Refine a triangle mesh every where using butterfly interpolation with coefficient w.

    *Example of use:*

    * `Refine all cells with butterfly interpolation (array) <Examples/Post/refine2.py>`_:

    .. literalinclude:: ../build/Examples/Post/refine2.py

    * `Refine all cells with butterfly interpolation (pyTree) <Examples/Post/refine2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/refine2PT.py

---------------------------------------

.. py:function:: Post.computeIndicatorValue (a, t, varName)

    Compute the indicator value on the unstructured octree mesh a based on the absolute maximum
    value of a varName field defined in the corresponding structured octree t.
    In the array version, t is a list of zones, and in the pyTree version, it can be a tree or a base or a list of bases
    or a zone or a list of zones.
    Variable varName can be located at nodes or centers.
    The resulting projected field is stored at centers in the octree mesh.

    *Example of use:*

    * `Project the maximum value of the indicator field on the octree mesh (array) <Examples/Post/.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeIndicatorValue.py

    * `Project the maximum value of the indicator field on the octree mesh (pyTree) <Examples/Post/PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeIndicatorValuePT.py

---------------------------------------

.. py:function:: Post.computeIndicatorField

    compute an indicator field to adapt an octree mesh with respect to the
    required number of points nbTargetPts, a field, and bodies.
    If refineFinestLevel=1, the finest level of the octree o is refined.
    If coarsenCoarsestLevel=1, the coarsest level of the octree o is
    coarsened provided the balancing is respected.<br>
    This function computes epsInf, epsSup, indicator such that when
    indicVal < valInf, the octree is coarsened (indicator=-1), when
    indicVal > valSup, the octree is refined (indicator=1).

    For an octree defined in an array o, and the field in indicVal:
    ::

        indicator, valInf, valSup = P.computeIndicatorField(o, indicVal, nbTargetPts=-1, bodies=[], refineFinestLevel=1, coarsenCoarsestLevel=1)

    For the pyTree version, the name varname of the field on which is based
    the indicator must be specified:
    ::

        o, valInf, valSup = P.computeIndicatorField(o, varname, nbTargetPts=-1, bodies=[], refineFinestLevel=1, coarsenCoarsestLevel=1)

    *Example of use:*

    * `Compute the adaptation indicator (array) <Examples/Post/computeIndicatorField.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeIndicatorField.py

    * `Compute the adaptation indicator (pyTree) <Examples/Post/computeIndicatorFieldPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeIndicatorFieldPT.py

---------------------------------------

Solution extraction
-------------------

.. py:function:: Post.extractPoint(A, (x,y,z), order=2, constraint=40., tol=1.e-6, hook=None)

    Extract the field in one or several points, given a solution defined by A.
    The extracted field(s) is returned as a list of values for each point.
    If the point (x,y,z) is not interpolable from a grid, then 0 for all fields is returned.

    To extract field in several points use:
    ::

        F = P.extractPoint(A, [(x1,y1,z1),(x2,y2,z2)], order=2, constraint=40., tol=1.e-6, hook=None)

    In the pyTree version, extractPoint returns the extracted solution
    from solutions located at nodes followed by the solution extracted from solutions at centers.

    If 'cellN', 'ichim', 'cellnf', 'status', or 'cellNF' variable is defined,
    it is returned in the last position in the output array.
    The interpolation order can be 2, 3, or 5.

    'constraint' is a thresold for extrapolation to occur. To enable more
    extrapolation, rise this value.

    If some blocks in A define surfaces, a tolerance 'tol' for interpolation cell search can be defined.

    A hook can be defined in order to keep in memory the ADT on the
    interpolation cell search. The hook argument must be a list of hooks, 
    each one being built for each donor zone using
    the "createHook" function of Converter module with 'extractMesh' argument.

    *Example of use:*

    * `Extraction in one point (array) <Examples/Post/extractPoint.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPoint.py

    * `Extraction in one point (pyTree) <Examples/Post/extractPointPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPointPT.py


---------------------------------------

.. py:function:: Post.extractPlane(A, (c1, c2, c3, c4), order=2, tol=1.e-6)

    slice a solution A with a plane.
    The extracted solution is interpolated from A.
    Interpolation order can be 2, 3, or 5
    (but the 5th order is very time-consuming for the moment).
    The best solution is kept. Plane is defined
    by :math:`c1\ x + c2\ y + c3\ z + c4 = 0`.

    *Example of use:*

    * `Extraction on a given plane (array) <Examples/Post/extractPlane.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPlane.py

    * `Extraction on a given plane (pyTree) <Examples/Post/extractPlanePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPlanePT.py


---------------------------------------

.. py:function:: Post.extractMesh(A, a, order=2, extrapOrder=1, constraint=40., tol=1.e-6, mode='robust', hook=None)

    Interpolate a solution from a set of donor zones defined by A to an extraction zone a.
    Parameter order can be 2, 3 or 5, meaning that 2nd, 3rd and 5th order interpolations are performed.

    Parameter 'constraint'>0 enables to extrapolate from A if interpolation is not possible for some points.
    Extrapolation order can be 0 or 1 and is defined by extrapOrder.

    If mode='robust', extract from the node mesh (solution in centers is first
    put to nodes, resulting interpolated solution is located in nodes).

    If mode='accurate', extract node solution from node mesh and center solution
    from center mesh (variables don't change location).

    The interpolation cell search can be preconditioned if extractMesh is applied several times using the same donor mesh.
    Parameter hook is only used in 'robust' mode and is a list of ADT (one per donor zone), each of them must be created and deleted by C.createHook and C.freeHook (see Converter module userguide).

    Exists also as in place version (_extractMesh) that modifies a and return None.

    *Example of use:*

    * `Extraction on an extraction zone (array) <Examples/Post/extractMesh.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractMesh.py

    * `Extraction on an extraction zone (pyTree) <Examples/Post/extractMeshPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractMeshPT.py

---------------------------------------

.. py:function:: Post.projectCloudSolution(pts, t, dim=3)

    Project the solution by a Least-Square Interpolation defined on a set of points pts defined as a 'NODE' zone
    to a body defined by a 'TRI' mesh in 3D and 'BAR' mesh in 2D.

    *Example of use:*

    * `projectCloudSolution (array) <Examples/Post/projectCloudSolution.py>`_:

    .. literalinclude:: ../build/Examples/Post/projectCloudSolution.py

    * `projectCloudSolution (pyTree) <Examples/Post/projectCloudSolutionPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/projectCloudSolutionPT.py


---------------------------------------

.. py:function:: Post.zipper(A, options=[])

    Build an unstructured unique surface mesh, given a list of structured
    overlapping surface grids A.
    Cell nature field is used to find blanked (0) and interpolated (2) cells.

    The options argument is a list of arguments such as ["argName", argValue]. Option names can be:

    - 'overlapTol' for tolerance required between two overlapping grids : if the projection distance between them is under this value then the grids are considered to be overset. Default value is 1.e-5.
    - For some cases, 'matchTol' can be set to modify the matching boundaries tolerance. Default value is set 1e-6.

    In most cases, one needn't modify this parameter.

    *Example of use:*

    * `Zipping of an overset surface (array) <Examples/Post/zipper.py>`_:

    .. literalinclude:: ../build/Examples/Post/zipper.py

    * `Zipping of an overset surface (pyTree) <Examples/Post/zipperPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/zipperPT.py


---------------------------------------

.. py:function:: Post.usurp(A)

    This function computes a "ratio" field for structured overlapping surfaces.
    The ratio field is located at cell centers. 
    In case of no overset, ratio are set to 1, otherwise ratio represents
    the percentage of overlap of a cell by another mesh.
    The finest cells have priority.
    All surfaces must be oriented in the same way.
    
    When using the array interface:
    ::

        C = P.usurp(A, B)

    the input arrays are a list of grid arrays A, defining nodes coordinates and a
    corresponding list of arrays defining the chimera nature of cells at cell centers B. Blanked cells must be flagged by a null value.
    Other values are equally considered as computed or interpolated cells.

    When using the pyTree interface:
    ::

        C = P.usurp(A)

    chimera cell nature field must be defined as a center field in A.

    Warning: normal of surfaces grids defined by A must be
    oriented in the same direction.

    *Example of use:*

    * `Ratio generation for the surface elements (array) <Examples/Post/usurp.py>`_:

    .. literalinclude:: ../build/Examples/Post/usurp.py

    * `Ratio generation for the surface elements (pyTree) <Examples/Post/usurpPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/usurpPT.py


---------------------------------------

Streams
-------

.. py:function:: Post.streamLine(A, (x0,y0,z0),  ['v1','v2,'v3'], N=2000, dir=2)

    Compute the stream line with N points starting from point (x0,y0,z0), given a solution A and a vector defined by 3 variables
    ['v1','v2,'v3'].
    Parameter 'dir' can be set to 1 (streamline follows velocity), -1
    (streamline follows -velocity), or 2
    (streamline expands in both directions).
    The output yields the set of N extracted points on the streamline,
    and the input fields at these points. The streamline computation
    stops when the current point is not interpolable from the input grids.

    *Example of use:*

    * `Streamline extraction (array) <Examples/Post/streamLine.py>`_:

    .. literalinclude:: ../build/Examples/Post/streamLine.py

    * `Streamline extraction (pyTree) <Examples/Post/streamLinePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/streamLinePT.py


---------------------------------------

.. py:function:: Post.streamRibbon(A, (x0,y0,z0), (nx,ny,nz), ['v1', 'v2', 'v3'], N=2000, dir=2)

0    Compute the stream ribbon starting from point (x0,y0,z0), of width and direction given by the vector (nx,ny,nz).
    This vector must be roughly orthogonal to the vector ['v1', 'v2', 'v3'] at point (x0,y0,z0).
    The output yields the set of N extracted points on the stream ribbon,
    and the input fields at these points. The stream ribbon computation
    stops when the current point is not interpolable from the input grids.

    *Example of use:*

    * `Stream ribbon extraction (array) <Examples/Post/streamRibbon.py>`_:

    .. literalinclude:: ../build/Examples/Post/streamRibbon.py

    * `Stream ribbon extraction (pyTree) <Examples/Post/streamRibbonPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/streamRibbonPT.py


---------------------------------------

.. py:function:: Post.streamSurf(A, c, ['v1','v2,'v3'], N=2000, dir=1)

    Compute the stream surface starting from a BAR array c.

    *Example of use:*

    * `Stream surface extraction (array) <Examples/Post/streamSurf.py>`_:

    .. literalinclude:: ../build/Examples/Post/streamSurf.py

    * `Stream surface extraction (pyTree) <Examples/Post/streamSurfPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/streamSurfPT.py


---------------------------------------

Isos
-------


.. py:function:: Post.isoLine(A, field, val)

    Compute an isoline correponding to value val of field.

    *Example of use:*

    * `Isoline computation (array) <Examples/Post/isoLine.py>`_:

    .. literalinclude:: ../build/Examples/Post/isoLine.py

    * `Isoline (pyTree) <Examples/Post/isoLinePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/isoLinePT.py


---------------------------------------

.. py:function:: Post.isoSurf(a, field, val, vars=None, split='simple')

    .. A1.O0.D0

    Compute an isosurface corresponding to value val of field (using marching
    tetrahedra). Resulting solution is always located in nodes.
    Return a list of two zones (one TRI and one BAR, if relevant).
    If vars (for ex: ['centers:F', 'G']) is given, extract only given variables.

    :param a:  input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param field: field name used in iso computation
    :type field: string
    :param val: value of field for extraction
    :type val: float
    :param vars: list of variable names you want to see on final iso-surface
    :type vars: list of strings
    :param split: 'simple' or 'withBarycenters', used in decomposing a in tetra (if needed)
    :type split: string

    *Example of use:*

    * `Isosurface extraction by marching tetra (array) <Examples/Post/isoSurf.py>`_:

    .. literalinclude:: ../build/Examples/Post/isoSurf.py

    * `Isosurface extraction by marching tetra (pyTree) <Examples/Post/isoSurfPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/isoSurfPT.py


---------------------------------------

.. py:function:: Post.isoSurfMC(a, field, val, vars=None, split='simple')

    .. A1.O0.D0

    Compute an isosurface corresponding to value val of field (using marching
    cubes). Resulting solution is always located in nodes.
    If vars (for ex: ['centers:F', 'G']) is given, extract only given variables.

    :param a: input data
    :type  a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param field: field name used in iso computation
    :type field: string
    :param val: value of field for extraction
    :type val: float
    :param vars: list of variable names you want to see on final iso-surface
    :type vars: list of strings
    :param split: 'simple' or 'withBarycenters', used in decomposing a in tetra (if needed)
    :type split: string
    :return: a list of isosurface (one per original zones)
    :rtype: list of arrays or list of zones

    
    *Example of use:*

    * `Isosurface by marching cube (array) <Examples/Post/isoSurfMC.py>`_:

    .. literalinclude:: ../build/Examples/Post/isoSurfMC.py

    * `Isosurface by marching cube (pyTree) <Examples/Post/isoSurfMCPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/isoSurfMCPT.py


---------------------------------------

Solution integration
--------------------

    For all integration functions, the interface is different when using
    Converter arrays interface or pyTree interface. For arrays, fields
    must be input separately, for pyTree, they must be defined in
    each zone.

---------------------------------------

.. py:function:: Post.integ(A, var='F')

    Compute the integral :math:`\int F.dS` of a scalar field (whose name is in var string) over
    the geometry defined by arrays containing the coordinates + field ( + an optional "ratio" field ).
    Solution and ratio can be located at nodes or at centers.

    For array interface:
    ::

        res = P.integ([coord], [field], [ratio]=[])

    For pyTree interface, the variable to be integrated can be specified. If no variable
    is specified, all the fields located at nodes and centers are integrated:
    ::

        res = P.integ(A, var='F')

    :param A: input data
    :type  A: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: field name
    :type var: string
    :return: the result of integration
    :rtype: a list of 1 float

    Exists also as parallel distributed version (P.Mpi.integ).

    *Example of use:*

    * `Scalar integration (array) <Examples/Post/integ.py>`_:

    .. literalinclude:: ../build/Examples/Post/integ.py

    * `Scalar integration (pyTree) <Examples/Post/integPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/integPT.py


---------------------------------------

.. py:function:: Post.integNorm(A, var='F')

    Compute the integral :math:`\int F.\vec n.dS` of a scalar field times the surface normal
    over the geometry defined by coord. For array interface:
    ::

        P.integNorm([coord], [field], [ratio]=[])

    For pyTree interface, the variable to be integrated can be specified. If no variable
    is specified, all the fields located at nodes and centers are integrated:
    ::

        P.integNorm(A, var='F')

    :param A: input data
    :type  A: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: field name
    :type var: string
    :return: the result of integration
    :rtype: double list of 3 floats

    Exists also as parallel distributed version (P.Mpi.integNorm).

    *Example of use:*

    * `Integration dot the surface normal (array) <Examples/Post/integNorm.py>`_:

    .. literalinclude:: ../build/Examples/Post/integNorm.py

    * `Integration dot the surface normal (pyTree) <Examples/Post/integNormPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/integNormPT.py


---------------------------------------

.. py:function:: Post.integNormProduct(A, vector=['vx','vy','vz'])

    Compute the integral :math:`\int \vec V \times \vec n.dS` of a vector field times the surface normal
    over the geometry defined by coord. The input field must have 3
    variables. For array interface, field must be a vector field:
    ::

        res = P.integNormProduct([coord], [field], [ratio]=[])

    For pyTree interface, the vector field to be integrated must be specified:
    ::

        res = P.integNormProduct(A, vector=['vx','vy','vz'])


    Exists also as parallel distributed version (P.Mpi.integNormProduct).

    *Example of use:*

    * `Integration cross the surface normal (array) <Examples/Post/integNormProduct.py>`_:

    .. literalinclude:: ../build/Examples/Post/integNormProduct.py

    * `Integration cross the surface normal (pyTree) <Examples/Post/integNormProductPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/integNormProductPT.py


---------------------------------------

.. py:function:: Post.integMoment(A, center=(0.,0.,0.), vector=['vx','vy','vz'])

    Compute the integral :math:`\int \vec{CM} \times \vec V.dS` of a moment over the geometry defined by coord.
    The input field must have 3 variables. center=(cx,cy,cz) are the center coordinates. 
    
    For array interface:
    ::

       res = P.integMoment([coord], [field], [ratio]=[], center=(0.,0.,0.))

    For pyTree interface, the vector of variables to be integrated must be specified:
    ::

       res = P.integMoment(A, center=(0.,0.,0.), vector=['vx','vy','vz'])

    Exists also as parallel distributed version (P.Mpi.integMoment).

    :param A: input data
    :type  A: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param vector: list of vector field names
    :type vector: list of 3 strings
    :return: the result of integration
    :rtype: a list of 3 floats

    *Example of use:*

    * `Moment integration  (array) <Examples/Post/integMoment.py>`_:

    .. literalinclude:: ../build/Examples/Post/integMoment.py

    * `Moment integration  (pyTree) <Examples/Post/integMomentPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/integMomentPT.py


---------------------------------------

.. py:function:: Post.integMomentNorm(A, center=(cx,cy,cz), var='F')

    Compute the integral :math:`\int \vec{CM} \times F.\vec n. dS` of a moment over the geometry 
    defined by coord, taking into account the surface normal. The input field is a scalar. 
    
    For array interface:
    ::

      res = P.integMomentNorm([coord], [field], [ratio]=[], center=(cx,cy,cz))

    For pyTree interface, the variable to be integrated can be specified. If no variable
    is specified, all the fields located at nodes and centers are integrated:
    ::

     res = P.integMomentNorm(A, center=(cx,cy,cz), var='F')

    :param A: input data
    :type  A: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: field name
    :type var: string
    :return: the result of integration
    :rtype: a list of 3 floats

    Exists also as parallel distributed version (P.Mpi.integMomentNorm).

    *Example of use:*

    * `Moment integration with normal (array) <Examples/Post/integMomentNorm.py>`_:

    .. literalinclude:: ../build/Examples/Post/integMomentNorm.py

    * `Moment integration with normal (pyTree) <Examples/Post/integMomentNormPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/integMomentNormPT.py

---------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



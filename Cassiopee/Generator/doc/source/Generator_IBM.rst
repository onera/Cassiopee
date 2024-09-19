.. Generator IBM documentation master file

:tocdepth: 2


Generator.IBM: mesh generation module for IBMs
==============================================================

Specific mesh generation functions for immersed boundaries (IB).

All of these functions can be executed in both sequential and parallel contexts.

.. py:module:: Generator.IBM


List of functions
#################


**-- IBM automatic grid generation**

.. autosummary::
   :nosignatures:

    Generator.IBM.buildOctree
    Generator.IBM.generateIBMMesh

   

Contents
#########
.. py:function:: Generator.IBM.generateIBMMesh(tb, dimPb=3, vmin=15, snears=0.01, dfars=10., tbox=None, to=None, octreeMode=0, check=False)

    Generates the full Cartesian mesh (octree/quadtree-based) for IBMs. The algorithm is divided into three main steps. It starts with the sequential octree generation from the surface definitions, through optional local adaptations from the refinement zones defined in tbox, to the resulting Cartesian mesh. The methodology is introduced and detailed in Peron and Benoit [2013, https://doi.org/10.1016/j.jcp.2012.07.029], and recalled in Constant [2023, http://dx.doi.org/10.13140/RG.2.2.35378.21449]. The resulting mesh is a collection of overset isotropic grids with minimal overlap.

    This function encapsulates the Generator.buildOctree function. For more details about the octree creation step, see the documentation of this function. If the octree has already been built, the user can also pass the octree as an input parameter (to). 

    This function fully operates in a distributed parallel environment and automatically splits the resulting Cartesian mesh into NP subzones, where NP is the number of MPI processes.

    :param tb: surface mesh
    :type tb: [zone, list of zones, base, tree]
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param vmin: minimum number of cells per direction for each octree level
    :type vmin: integer
    :param snears: minimum cell spacing(s) near the bodies
    :type snears: float or list of floats
    :param dfars: extent(s) of the domain from the bodies
    :type dfars: float or list of floats
    :param tbox: refinement bodies
    :type tbox: [zone, list of zones, base, tree]
    :param to: input octree if already created
    :type to: [zone, list of zones, base, tree]
    :param octreeMode: octree generation mode
    :type octreeMode: 0 or 1
    :param check: if True: write octree.cgns locally
    :type check: boolean
    :return: block-structured mesh tree

    *Example of use:*

    * `Generates the full Cartesian mesh for IBMs (pyTree) <Examples/Generator/generateIBMMeshPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/generateIBMMeshPT.py

---------------------------------------

.. py:function:: Generator.IBM.buildOctree(tb, dimPb=3, vmin=15, snears=0.01, dfars=10., tbox=None, octreeMode=0)

    Builds an octree (3D) or quadtree (2D) tree from the surface definitions stored in tb. This function is inherently sequential, and the geometry file must be shared among all processors when running in parallel. The resulting octree (or quadtree) is balanced to respect a maximum ratio of 2 between adjacent leaf nodes. By default, the current balancing mode also respects the same condition on nodes connected by one vertice. 

    Since this function, which is based on the Generator.PyTree.octree() function, is primarily used to automatically generate Cartesian grids around immersed boundaries, a final expansion of the lowest level leaf nodes is performed so that the minimum spacing imposed near the wall is sufficiently propagated in the wall normal direction. Local refinement zones stored in the tbox argument can be used to further refine the octree.

    This function takes into account three main parameters which are vmin, snears and dfars. 
    
    * vmin is a global parameter that controls the minimum number of cells per direction for each octree level. For example, a small vmin value will result in a small number of points against a large number of elementary Cartesian blocks, as well as more frequent resolution changes from the wall boundaries. Empirically, one should use vmin values between 7 and 21.

    * snears defines the minimum near-wall spacing. This argument can be passed globally as a float, or locally as a list of floats whose size must be equal to the number of zones in the tb file. Note that these values will eventually be overwritten by any snears values found in each subzone of the geometry pytree (see Geom.IBM.setSnear).

    * The dfars argument specifies the global domain extent from the geometry bounding boxes. Like snears, this argument can be a float or a list of floats, and local values found in tb will eventually overwrite the values passed as argument (see Geom.IBM.setDfar).

    Since the octree is created by recursively subdividing cubes into octants, only the final snear or dfar values can be exact. The parameter octreeMode allows the user to generate an octree by fixing one or the other. By default, octreeMode is set to 0 and the domain extent is fixed. The subdivision step ends when the minimum near-wall spacing is close enough to the minimum snear value specified by the user. Note that  in some cases the actual snear can be up to 20% lower or higher than the expected snear value(s). When octreeMode is set to 1, the minimum near-wall spacing is fixed and the domain extent is finally modified to get as close as possible to the desired dfars values.

    :param tb: surface mesh
    :type tb: [zone, list of zones, base, tree]
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param vmin: minimum number of cells per direction for each octree level
    :type vmin: integer
    :param snears: minimum cell spacing(s) near the bodies
    :type snears: float or list of floats
    :param dfars: extent(s) of the domain from the bodies
    :type dfars: float or list of floats
    :param tbox: refinement bodies
    :type tbox: [zone, list of zones, base, tree]
    :param octreeMode: octree generation mode
    :type octreeMode: 0 or 1
    :return: monozone octree (3D) or quadtree (2D), Quad (2D) or Hex (3D) type

    *Example of use:*

    * `Builds an octree from the surface definitions (pyTree) <Examples/Generator/buildOctreePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/buildOctreePT.py
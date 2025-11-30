.. Geom documentation master file

:tocdepth: 2

Geom: Geometry Definition Module
=========================================


Preamble
########

In this module, a geometry is defined discretely with a great number of points.
The discrete geometry is stored in a Converter array (as defined in Converter documentation) and in
a CGNS/Python tree zone or set of zones, depending on the selected interface. 

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

To use the module with the Converter array interface::

   import Geom as D

To use the module with the CGNS/Python interface::

   import Geom.PyTree as D


.. py:module:: Geom

List of functions
#################

**-- Geometry creation**

.. autosummary::
   :nosignatures:

   Geom.point
   Geom.cloud
   Geom.line
   Geom.polyline
   Geom.circle
   Geom.naca
   Geom.profile
   Geom.spline
   Geom.nurbs
   Geom.bezier
   Geom.curve
   Geom.surface
   Geom.cone
   Geom.torus
   Geom.sphere
   Geom.sphere6
   Geom.sphereYinYang
   Geom.disc
   Geom.triangle
   Geom.quadrangle
   Geom.box
   Geom.cylinder

**-- Typing text using meshes**

.. autosummary::
   :nosignatures:

   Geom.text1D
   Geom.text2D
   Geom.text3D

**-- Geometry modification**

.. autosummary::
   :nosignatures:

   Geom.uniformize
   Geom.refine
   Geom.enforceh 
   Geom.lineDrive
   Geom.orthoDrive
   Geom.axisym
   Geom.connect1D
   Geom.Offset.offsetSurface

   .. Geom.volumeFromCrossSections
   .. Geom.addSeparationLine

**-- 1D distributions**

.. autosummary::
   :nosignatures:

    Geom.distrib1
    Geom.distrib2
    Geom.getDistribution
   
**-- Information about geometries**

.. autosummary::
   :nosignatures:

   Geom.getLength
   Geom.getDistantIndex
   Geom.getNearestPointIndex
   Geom.getCurvatureRadius
   Geom.getCurvatureAngle
   Geom.getCurvatureHeight
   Geom.getSharpestAngle
   Geom.getCurvilinearAbscissa
   Geom.getTangent
   Geom.getUV
   Geom.getUVFromIJ

Contents
########

Geometry creation
-----------------
A geometry can be defined by either structured (i-arrays in 1D, (i,j) arrays in 2D) or unstructured 
(BAR arrays in 1D and QUAD or TRI arrays in 2D) grids.
 
A polyline is defined as a C0 i-array which contains only the polyline points (with no extra discretization points).

---------------------------------------

.. py:function:: Geom.point(P)

    Create a point of coordinates P=(x,y,z).

    :param P: (x,y,z) of point
    :type P: 3-tuple of floats 
    :return: a point
    :rtype: one array/zone (NODE)

    *Example of use:*

    * `Creation of a point (array) <Examples/Geom/point.py>`_:

    .. literalinclude:: ../build/Examples/Geom/point.py

    * `Creation of a point (pyTree) <Examples/Geom/pointPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/pointPT.py
    
---------------------------------------

.. py:function:: Geom.cloud(arr)

    Create a point cloud of coordinates arr=(X,Y,Z) where X=[x1, x2, ..., xn]

    :param arr: (x,y,z) of point cloud
    :type arr: 3-tuple of lists/numpy.ndarrays
    :return: a point cloud
    :rtype: one array/zone (NODE)

    *Example of use:*

    * `Creation of a point cloud (array) <Examples/Geom/cloud.py>`_:

    .. literalinclude:: ../build/Examples/Geom/cloud.py

    * `Creation of a point cloud (pyTree) <Examples/Geom/cloudPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/cloudPT.py

---------------------------------------

.. py:function:: Geom.line(P1,P2,N=100)

    Create a line from point P1 to point P2, uniformly discretized with N points.

    :param P1:  (x,y,z) of the starting point
    :type P1: 3-tuple of floats 
    :param P2:  (x,y,z) of the end point
    :type P2: 3-tuple of floats 
    :param N:  number of points discretizing the line
    :type N: integer 
    :return: a line
    :rtype: one array/zone (1D STRUCT)

    *Example of use:*

    * `Creation of a line (array) <Examples/Geom/line.py>`_:

    .. literalinclude:: ../build/Examples/Geom/line.py

    * `Creation of a line (pyTree) <Examples/Geom/linePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/linePT.py

---------------------------------------

.. py:function:: Geom.polyline(Pts)

    Create a polyline made of points Pts=[P1, P2,...,Pn].

    :param Pts:  list of (x,y,z) of points defining the polyline
    :type Pts: list of 3-tuple of floats 
    :return: a polyline
    :rtype: one array/zone (1D STRUCT)

    *Example of use:*

    * `Creation of a polyline (array) <Examples/Geom/polyline.py>`_:

    .. literalinclude:: ../build/Examples/Geom/polyline.py

    * `Creation of a polyline (pyTree) <Examples/Geom/polylinePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/polylinePT.py

---------------------------------------

.. py:function:: Geom.circle(C, R, tetas=0., tetae=360., N=100)

    Create a circle or a circle arc of center C and radius R.

    :param C:  (x,y,z) of circle center
    :type C: 3-tuple of floats 
    :param R:  radius of the circle
    :type R: float 
    :param tetas:  initial azimuth of the circle arc
    :type tetas: float
    :param tetae: end azimuth of circle arc
    :type tetae: float
    :param N: number of points discretizing the circle arc    
    :type N: integer
    :return: a circle
    :rtype: one array/zone (1D STRUCT)

    *Example of use:*

    * `Creation of a circle (array) <Examples/Geom/circle.py>`_:

    .. literalinclude:: ../build/Examples/Geom/circle.py

    * `Creation of a circle (pyTree) <Examples/Geom/circlePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/circlePT.py

---------------------------------------

.. py:function:: Geom.naca(e, N=101, sharpte=True)

    Create a NACA profile. e can be the thickness of the profile 
    (e=15. for NACA0015 for instance) or a string of digit such as
    "0012" for serie 4, "23012" for serie 5, "0008-45" for modified
    serie 4 of NACA profiles.

    :param e:  thickness of the NACA00xx profile or digit string
    :type e: float or string
    :param N: number of points discretizing the profile    
    :type N: integer
    :param sharpte: true if sharp trailing edge
    :type sharpte: boolean
    :return: a NACAxx profile
    :rtype: one array/zone (1D STRUCT)

    *Example of use:*

    * `Creation of a NACA0012 (array) <Examples/Geom/naca.py>`_:

    .. literalinclude:: ../build/Examples/Geom/naca.py

    * `Creation of a NACA0012 (pyTree) <Examples/Geom/nacaPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/nacaPT.py

---------------------------------------

.. py:function:: Geom.profile(name=None)

    Create a wing profile of given name. 
    The name must be of type "series/profile".
    If name is not provided, list available profiles.

    :param name: profile name
    :type name: string of type "series/profile"
    :return: a profile
    :rtype: one array/zone (1D STRUCT)

    *Example of use:*

    * `Creation of a profile (array) <Examples/Geom/profile.py>`_:

    .. literalinclude:: ../build/Examples/Geom/profile.py

    * `Creation of a profile (pyTree) <Examples/Geom/profilePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/profilePT.py

---------------------------------------

.. py:function:: Geom.spline(Pts, order=3, N=100, M=100, density=-1)

    Create a Spline curve/surface using control points defined by Pts.

    :param Pts: i-mesh (resp. (i,j)-mesh) of control points for a Spline curve (resp. surface)
    :type Pts: array or zone of control points
    :param order: order of the Spline
    :type order: integer
    :param N: number of points in the i-direction in the resulting discretized Spline
    :type N: integer
    :param M: number of points in the j-direction in the resulting discretized Spline
    :type M: integer
    :param density: density of points in the discretized Spline  (instead of specifying N and M)  
    :type density: float
    :return: a Spline curve/surface
    :rtype: one array/zone (1D STRUCT or 2D STRUCT)

    *Example of use:*

    * `Creation of a Spline (array) <Examples/Geom/spline.py>`_:

    .. literalinclude:: ../build/Examples/Geom/spline.py

    * `Creation of a Spline (pyTree) <Examples/Geom/splinePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/splinePT.py

---------------------------------------

.. py:function:: Geom.nurbs(Pts, W, order=3, N=100, M=100, density=-1.)

    Create a NURBS curve/surface using control points and weights defined by Pts and W.

    :param Pts: i-mesh (resp. (i,j)-mesh) of control points for a NRUBS curve (resp. surface)
    :type Pts: array or zone
    :param W: weight for each control point defined in Pts
    :type W: string 
    :param order: order of the NURBS    
    :type order: integer
    :param N: number of points in the i-direction in the resulting discretized NURBS
    :type N: integer
    :param M: number of points in the j-direction in the resulting discretized NURBS
    :type M: integer
    :param density: density of points in the discretized NURBS (instead of specifying N and M)  
    :type density: float
    :return: a NURBS curve/surface
    :rtype: one array/zone (1D STRUCT or 2D STRUCT)

    *Example of use:*

    * `Creation of a NURBS (array) <Examples/Geom/nurbs.py>`_:

    .. literalinclude:: ../build/Examples/Geom/nurbs.py

    * `Creation of a NURBS (pyTree) <Examples/Geom/nurbsPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/nurbsPT.py

---------------------------------------

.. py:function:: Geom.bezier(Pts, N=100, M=100, density=-1.)

    Create a Bezier curve/surface using control points defined by Pts.

    :param Pts: i-mesh (resp. (i,j)-mesh) of control points for a spline curve (resp. surface)
    :type Pts: array or zone
    :param N: number of points in the i-direction in the resulting discretized Bezier
    :type N: integer
    :param M: number of points in the j-direction in the resulting discretized Bezier
    :type M: integer
    :param density: density of points in the discretized Bezier (instead of specifying N and M)  
    :type density: float
    :return: a Bezier curve/surface
    :rtype: one array/zone (1D STRUCT or 2D STRUCT)

    *Example of use:*

    * `Creation of a Bezier curve (array) <Examples/Geom/bezier.py>`_:

    .. literalinclude:: ../build/Examples/Geom/bezier.py

    * `Creation of a Bezier curve (pyTree) <Examples/Geom/bezierPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/bezierPT.py

---------------------------------------

.. py:function:: Geom.curve(f, N=100)

    Create a curve defined by a parametric function or an expression.

    :param f: Python function or set of expressions separated by ";""
    :type f: Python function or string
    :param N: number of discretization points per direction
    :type N: integer
    :return: a parametric curve
    :rtype: one array/zone (1D STRUCT)

    *Example of use:*

    * `Creation of a parametric curve (array) <Examples/Geom/curve.py>`_:

    .. literalinclude:: ../build/Examples/Geom/curve.py

    * `Creation of a parametric curve (pyTree) <Examples/Geom/curvePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/curvePT.py

---------------------------------------

.. py:function:: Geom.surface(f, N=100)

    Create a surface defined by an parametric function or an expression.

    :param f: Python function or set of expressions separated by ";""
    :type f: Python function or string
    :param N: number of discretization points per direction
    :type N: integer
    :return: a parametric surface
    :rtype: one array/zone (2D STRUCT)

    *Example of use:*

    * `Creation of a parametric surface (array) <Examples/Geom/surface.py>`_:

    .. literalinclude:: ../build/Examples/Geom/surface.py

    * `Creation of a parametric surface (pyTree) <Examples/Geom/surfacePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/surfacePT.py

---------------------------------------

.. py:function:: Geom.cone(C, Rb, Rt, H, N=100)

    Create a cone discretized by NxN points.

    :param C: center coordinates
    :type C: 3-tuple of floats
    :param Rb: radius of the basis of the cone
    :type Rb: float
    :param Rt: radius of the top of the cone
    :type Rt: float
    :param H: height of the cone
    :type H: float
    :param N: number of discretization points per direction
    :type N: integer
    :return: a cone
    :rtype: one array/zone (2D STRUCT)

    *Example of use:*

    * `Creation of a cone(array) <Examples/Geom/cone.py>`_:

    .. literalinclude:: ../build/Examples/Geom/cone.py

    * `Creation of a cone (pyTree) <Examples/Geom/conePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/conePT.py

---------------------------------------

.. py:function:: Geom.torus(C,R,r,alphas=0.,alphae=360.,betas=0.,betae=360.,NR=100,Nr=100)

    Create a portion of torus discretized by NRxNr points, of center C, axis Z and radii R (main radius) and r (tube radius) between angles
    alphas and alphae (in the XY-plane) and between betas and betae (in the RZ-plane).

    :param C: center coordinates
    :type C: 3-tuple of floats
    :param R: main radius
    :type R: float
    :param r: tube radius
    :type r: float
    :param alphas: minimum azimuth in the XY-plane
    :type alphas: float
    :param alphae: maximum azimuth in the XY-plane
    :type alphae: float    
    :param betas: minimum azimuth in the RZ-plane
    :type betas: float
    :param betae: maximum azimuth in the RZ-plane
    :type betae: float    
    :param NR: number of discretization points in azimuth
    :type NR: integer
    :param Nr: number of discretization points in the axial direction
    :type Nr: integer
    :return: a torus
    :rtype: one array/zone (2D STRUCT)

    *Example of use:*

    * `Creation of a torus (array) <Examples/Geom/torus.py>`_:

    .. literalinclude:: ../build/Examples/Geom/torus.py

    * `Creation of a torus (pyTree) <Examples/Geom/torusPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/torusPT.py

---------------------------------------

.. py:function:: Geom.sphere(C, R, N=100)

    Create a structured mesh defining a sphere of radius R with N points in the longitudinal direction and 2xN along the latitude.

    :param C: sphere center coordinates
    :type C: 3-tuple of floats
    :param R: sphere radius
    :type R: float
    :param N: number of discretization points in the longitudinal direction
    :type N: integer
    :return: a structured mesh of a sphere degenerated at poles
    :rtype: one array/zone (2D STRUCT)

    *Example of use:*

    * `Creation of a sphere (array) <Examples/Geom/sphere.py>`_:

    .. literalinclude:: ../build/Examples/Geom/sphere.py

    * `Creation of a sphere (pyTree) <Examples/Geom/spherePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/spherePT.py

---------------------------------------

.. py:function:: Geom.sphere6(C, R, N=100, ntype='STRUCT')

    Create a mesh made of 6 parts defining a sphere of radius R with N points per direction.
    This mesh is not degenerated at poles in consequence.

    :param C: sphere center coordinates
    :type C: 3-tuple of floats
    :param R: sphere radius
    :type R: float
    :param N: number of discretization points in the longitudinal direction
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a sphere
    :rtype: a list of 6 arrays/zones (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Non-degenerated surface mesh of a sphere (array) <Examples/Geom/sphere6.py>`_:

    .. literalinclude:: ../build/Examples/Geom/sphere6.py

    * `Non-degenerated surface mesh of a sphere (pyTree) <Examples/Geom/sphere6PT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/sphere6PT.py

---------------------------------------

.. py:function:: Geom.sphereYinYang(C, R, N=100, ntype='STRUCT')

    Create an overset mesh of 2 parts defining a sphere of radius R with N points per direction.

    :param C: sphere center coordinates
    :type C: 3-tuple of floats
    :param R: sphere radius
    :type R: float
    :param N: number of discretization points in the longitudinal direction
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a sphere
    :rtype: a list of 2 arrays/zones (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Creation of a Yin-Yang sphere (array) <Examples/Geom/sphereYinYang.py>`_:

    .. literalinclude:: ../build/Examples/Geom/sphereYinYang.py

    * `Creation of a Yin-Yang sphere (pyTree) <Examples/Geom/sphereYinYangPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/sphereYinYangPT.py

---------------------------------------

.. py:function:: Geom.disc(C, R, N=100, ntype='STRUCT')

    Create a mesh of 5 parts defining a disc of radius R with NxN grid points.

    :param C: sphere center coordinates
    :type C: 3-tuple of floats
    :param R: disc radius
    :type R: float
    :param N: number of discretization points for each grid
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a disc
    :rtype: a list of 5 arrays/zones (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Creation of a disc (array) <Examples/Geom/disc.py>`_:

    .. literalinclude:: ../build/Examples/Geom/disc.py

    * `Creation of a disc (pyTree) <Examples/Geom/discPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/discPT.py

---------------------------------------

.. py:function:: Geom.triangle(P1, P2, P3, N=0, ntype='TRI')

    Create a triangle mesh defined by 3 vertices P1, P2, P3.

    :param P1:  (x,y,z) of first vertex
    :type P1: 3-tuple of floats 
    :param P2:  (x,y,z) of second vertex
    :type P2: 3-tuple of floats 
    :param P3:  (x,y,z) of third vertex
    :type P3: 3-tuple of floats
    :param N: number of discretization points
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a triangle
    :rtype: a list of 3 arrays/zones (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Creation of a triangle (array) <Examples/Geom/triangle.py>`_:

    .. literalinclude:: ../build/Examples/Geom/triangle.py

    * `Creation of a triangle (pyTree) <Examples/Geom/trianglePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/trianglePT.py

---------------------------------------

.. py:function:: Geom.quadrangle(P1, P2, P3, P4, N=0, ntype='QUAD')

    Create a quadrangle of vertices P1, P2, P3, P4.

    :param P1:  (x,y,z) of first vertex
    :type P1: 3-tuple of floats 
    :param P2:  (x,y,z) of second vertex
    :type P2: 3-tuple of floats 
    :param P3:  (x,y,z) of third vertex
    :type P3: 3-tuple of floats     
    :param P4:  (x,y,z) of fourth vertex
    :type P4: 3-tuple of floats
    :param N: number of discretization points
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a quadrangle
    :rtype: a list of 1 array/zone (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Creation of a quadrangle (array) <Examples/Geom/quadrangle.py>`_:

    .. literalinclude:: ../build/Examples/Geom/quadrangle.py

    * `Creation of a quadrangle (pyTree) <Examples/Geom/quadranglePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/quadranglePT.py


-------------------------------------------------------------------

.. py:function:: Geom.box(P1, P2, N=100, ntype='STRUCT')

    Create an axis aligned box passing by points P1 and P2.

    :param P1:  (x,y,z) of first vertex
    :type P1: 3-tuple of floats 
    :param P2:  (x,y,z) of second vertex
    :type P2: 3-tuple of floats 
    :param N: number of discretization points
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a box
    :rtype: a list of 6 arrays/zones (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Creation of a box (array) <Examples/Geom/box.py>`_:

    .. literalinclude:: ../build/Examples/Geom/box.py

    * `Creation of a box (pyTree) <Examples/Geom/boxPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/boxPT.py

-------------------------------------------------------------------

.. py:function:: Geom.cylinder(C, R, H, N=100, ntype='STRUCT')

    Create cylinder mesh made of two discs of center C and radius R and
    of height H.

    :param C:  bottom disc center
    :type C: 3-tuple of floats 
    :param R:  Radius of discs
    :type R: float
    :param H: Height of cylinder
    :type H: float
    :param N: number of discretization points
    :type N: integer
    :param ntype: type of output mesh ('STRUCT', 'QUAD', 'TRI')
    :type ntype: string
    :return: a mesh of a cylinder
    :rtype: a list of 11 arrays/zones (STRUCT) or 1 array/zone (QUAD and TRI)

    *Example of use:*

    * `Creation of a cylinder (array) <Examples/Geom/cylinder.py>`_:

    .. literalinclude:: ../build/Examples/Geom/cylinder.py

    * `Creation of a cylinder (pyTree) <Examples/Geom/cylinderPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/cylinderPT.py

---------------------------------------------------------------------------


Typing text using meshes
------------------------

.. py:function:: Geom.text1D(text, font='vera', smooth=0, offset=0.5)

    Create 1D meshes of given text.

    :param text:   text with separated characters
    :type text:    string
    :param font:   chosen font name (can be 'vera','chancery','courier','text1','nimbus')
    :type font:    string
    :param smooth: letter smoothness (0-4)
    :type smooth:  integer
    :param offset: distance between two letters
    :type offset:  float
    :return:       a mesh for each character of the text
    :rtype:        a list of arrays or zones

    *Example of use:*

    * `Text defined by a set of 1D meshes (array) <Examples/Geom/text1D.py>`_:

    .. literalinclude:: ../build/Examples/Geom/text1D.py

    * `Text defined by a set of 1D meshes (pyTree) <Examples/Geom/text1DPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/text1DPT.py

---------------------------------------

.. py:function:: Geom.text2D(text, font='vera', smooth=0, offset=0.5)

    Create a triangular mesh of a text (letters are filled with triangles).

    :param text: text with separated characters
    :type text: string
    :param font: chosen font name (can be 'vera','chancery','courier','text1','nimbus')
    :type font: string
    :param smooth: letter smoothness (0-4)
    :type smooth: integer
    :param offset: distance between two letters
    :type offset: float
    :return: a single mesh for the text string
    :rtype: an array or a zone

    *Example of use:*

    * `Text defined by a set of 2D meshes (array) <Examples/Geom/text2D.py>`_:

    .. literalinclude:: ../build/Examples/Geom/text2D.py

    * `Text defined by a set of 2D meshes (pyTree) <Examples/Geom/text2DPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/text2DPT.py

---------------------------------------

.. py:function:: Geom.text3D(text,font='vera',smooth=0,offset=0.5,thickness=8.)

    Create a 3D mesh of a text.

    :param text: text with separated characters
    :type text: string
    :param font: chosen font name (can be 'vera','chancery','courier','text1','nimbus')
    :type font: string
    :param smooth: letter smoothness (0-4)
    :type smooth: integer
    :param offset: distance between two letters
    :type offset: float
    :param thickness: thickness of letters
    :type thickness: float
    :return: a single mesh of text
    :rtype: an array or a zone

    *Example of use:*

    * `Text defined by a set of 3D meshes (array) <Examples/Geom/text3D.py>`_:

    .. literalinclude:: ../build/Examples/Geom/text3D.py

    * `Text defined by a set of 3D meshes (pyTree) <Examples/Geom/text3DPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/text3DPT.py

---------------------------------------------------------------------------

Geometry modification
------------------------

.. py:function:: Geom.uniformize(a, N=100, h=-1, factor=-1, density=-1, sharpAngle=30.)

    Remesh a 1D curve with a regular mesh step.
    You can specify one of N or factor or density or h.

    :param a: original curve to be remeshed (i-STRUCT or BAR) 
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param N: the final number of points
    :type N: int
    :param h: the final mesh step
    :type h: float
    :param factor: factor for the number of points regarding initial number of points of curve.
    :type factor: float
    :param density: point density
    :type density: float
    :param sharpAngle: point where the curve has a local angle greated than sharpAngle are enforced.
    :type sharpAngle: float
    :rtype: identical to a

    *Example of use:*

    * `Uniformizing steps of a 1D curve (array) <Examples/Geom/uniformize.py>`_:

    .. literalinclude:: ../build/Examples/Geom/uniformize.py

    * `Uniformizing steps of a 1D curve (pyTree) <Examples/Geom/uniformizePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/uniformizePT.py

---------------------------------------------------------------

.. py:function:: Geom.refine(a, N=10, factor=-1, sharpAngle=30.)

    Remesh a 1D curve keeping the original point distribution but densifying or
    coarsening it. You can specify N or factor.

    :param a: original curve to be refined (i-STRUCT or BAR)
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param N: the final number of points
    :type N: int
    :param factor: factor for the number of points regarding initial number of points of curve.
    :type factor: float
    :param sharpAngle: point where the curve has a local angle greated than sharpAngle are enforced.
    :type sharpAngle: float
    :rtype: identical to a

    *Example of use:*

    * `Refining/coarsening a 1D curve (array) <Examples/Geom/refine.py>`_:

    .. literalinclude:: ../build/Examples/Geom/refine.py

    * `Refining/coarsening a 1D curve (pyTree) <Examples/Geom/refinePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/refinePT.py

---------------------------------------------------------------

.. py:function:: Geom.enforceh(a, N=100, h=-1)

    Enforce some mesh steps or mesh factors in a 1D curve.
    To enforce a step use D.setH(a, ind, h), to enforce a factor, 
    use D.setF(a, ind, f).
    If you want to enforce h, you must specify N, the final number of points.
    If you want to enforce f, you must specify h, the mesh size for f=1.

    :param a: original curve to be remeshed (i-STRUCT or BAR)
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param N: the final number of points
    :type N: int
    :param h: the final mesh step
    :type h: float
    :rtype: identical to a

    *Example of use:*

    * `Enforce steps in a 1D curve (array) <Examples/Geom/enforceh.py>`_:

    .. literalinclude:: ../build/Examples/Geom/enforceh.py

    * `Enforce steps in a 1D curve (pyTree) <Examples/Geom/enforcehPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/enforcehPT.py

---------------------------------------------------------------

.. py:function:: Geom.lineDrive(a, d)

    Generate a surface mesh starting from a curve a and a single or a set of driving curves.
    The first point of the driving curves must match with one point of the original curve a.

    :param a: original curve to be extruded wrt the driving curve d 
    :type a: array or zone 
    :param d: driving curve or set of driving curves
    :type d: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :return: a surface structured mesh
    :rtype: an array or a zone

    *Example of use:*

    * `Extrusion with driving curve of a NACA0012 (array) <Examples/Geom/lineDrive.py>`_:

    .. literalinclude:: ../build/Examples/Geom/lineDrive.py

    * `Extrusion with driving curve of a NACA0012 (pyTree) <Examples/Geom/lineDrivePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/lineDrivePT.py

---------------------------------------------------------------

.. py:function:: Geom.orthoDrive(a, d, mode=0)

    Generate a surface mesh starting from a curve a and a single driving curve.
    The initial mesh is driven orthogonally to the driving curve.
    The first point of the driving curves must match with one point of the original curve a.

    :param a: original curve to be extruded orthogonally wrt the driving curve d 
    :type a: array or zone 
    :param d: driving curve
    :type d: [array] or [zone]
    :param mode: if mode=0, return one single zone, if mode=1, duplicate a and return a list of zones
    :return: a surface structured mesh or a list of meshes

    *Example of use:*

    * `Orthogonal extrusion of a NACA0012 (array) <Examples/Geom/orthoDrive.py>`_:

    .. literalinclude:: ../build/Examples/Geom/orthoDrive.py

    * `Orthogonal extrusion of a NACA0012 (pyTree) <Examples/Geom/orthoDrivePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/orthoDrivePT.py

.. ---------------------------------------------------------------
.. .. py:function:: Geom.addSeparationLine(a, b)
    Add a separation line between two geometries in contact.
    :param a: first geometry
    :type a: array or zone 
    :param b: separation line
    :type b: array or zone
    :return: two geometries with b defined in both geometries
    :rtype: list of arrays or list of zones
    *Example of use:*
    * `Separation line after a circle (array) <Examples/Geom/addSeparationLine.py>`_:
    .. literalinclude:: ../build/Examples/Geom/addSeparationLine.py
    * `Separation line after a circle (pyTree) <Examples/Geom/addSeparationLinePT.py>`_:
    .. literalinclude:: ../build/Examples/Geom/addSeparationLinePT.py

---------------------------------------

.. py:function:: Geom.axisym(a, C, axis, angle=360., Ntheta=100, rmod=None)

    Create an axisymmetrical mesh given one of its borders following axis.
    
    Exists also as in place version (_axisym) that modifies a and returns None.

    :param a: axis-aligned border of the axisymmetrical mesh (either structured-1D or 2D or BAR or TRI or QUAD)
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param C: center of rotation of the mesh
    :type C: 3-tuple of floats
    :param axis: rotation axis
    :type axis: 3-tuple of floats
    :param angle: azimuthal sector angle
    :type angle: float
    :param Ntheta: number of points in the azimuthal direction
    :type Ntheta: integer
    :param rmod: optional curve defining r=f(theta) instead of defining theta and Ntheta
    :type rmod: identical to a
    :return: a 2D or 3D mesh (either structured or QUAD or PENTA or HEXA)
    :rtype: Identical to a

    *Example of use:*

    * `Creation of a surface mesh by axisymmetry  (array) <Examples/Geom/axisym.py>`_:

    .. literalinclude:: ../build/Examples/Geom/axisym.py

    * `Creation of a surface mesh by axisymmetry (pyTree) <Examples/Geom/axisymPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/axisymPT.py


.. ---------------------------------------------------------------
.. .. py:function:: Geom.volumeFromCrossSections(contours)
    Create a volume mesh from cross-sections curves in the (X,Y) plane. 
    One curve is required per cross-section at different heights Z.
    :param contours: set of cross-section curves
    :type contours: list of arrays or list of zones 
    :return: the volume of the solid body
    :rtype: array or zone
    *Example of use:*
    * `Creation of a volume mesh from contours (array) <Examples/Geom/volumeFromCrossSection.py>`_:
    .. literalinclude:: ../build/Examples/Geom/volumeFromCrossSections.py
    * `Creation of volume mesh from contours (pyTree) <Examples/Geom/volumeFromCrossSectionsPT.py>`_:
    .. literalinclude:: ../build/Examples/Geom/volumeFromCrossSectionsPT.py

---------------------------------------

.. py:function:: Geom.connect1D(curves,sharpness=0,N=10,lengthFactor=1.)
  
    Connect non-matching curves by a line or by a Spline with N points.

    :param curves: two curves to be connected
    :type curves: list of arrays or list of zones 
    :param sharpness: 0: connected by a line; 1: connected by a Spline
    :type sharpness: integer
    :param N: number of points in the connection
    :type N: integer
    :param lengthFactor: the connection is bounded by lengthFactor x the length of the initial curves. 
    :type lengthFactor: float
    :return: a single curve connecting both curves
    :rtype: array or zone

    *Example of use:*

    * `Connect two lines  (array) <Examples/Geom/connect1D.py>`_:

    .. literalinclude:: ../build/Examples/Geom/connect1D.py

    * `Connect two lines (pyTree) <Examples/Geom/connect1DPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/connect1DPT.py

---------------------------------------

.. py:function:: Geom.Offset.offsetSurface(a, offset=1., pointsPerUnitLength=1., algo=0, dim=3)
  
    Offset a surface of certain distance. The underlaying algorithm uses an 
    octree (algo=1) or a cartesian grid (algo=0) with a number of points normalized per unit length
    determined by the parameter pointsPerUnitLength. Distance is then computed on this grid 
    and an isosurface of given "offset" distance is output.

    :param a: surface to be offset
    :type a: list of arrays or list of zones 
    :param offset: distance where the offset surface is created
    :type offset: float
    :param pointsPerUnitLength: number of points  per unit length
    :type pointsPerUnitLength: float
    :param algo: 0: cartesian (faster); 1: octree (less memory)
    :type algo: integer 
    :param dim: dimension of the pb (2 or 3)
    :type dim: integer 
    :return: a single zone
    :rtype: array or zone

    *Example of use:*


    * `Offset a surface (pyTree) <Examples/Geom/offsetSurfacePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/offsetSurfacePT.py


---------------------------------------------------------------------------

1D distributions
------------------

.. py:function:: Geom.distrib1(a, h, normalized=True)

    Return the distribution enforcing h everywhere in curve a.

    :param a: 1D mesh
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones] 
    :param h: size to enforce
    :type h: float
    :param normalized: if True, return a normalized distribution
    :type normalized: boolean
    :return: the distribution of the curve as mesh coordinates with h enforced
    :rtype: Identical to a

    *Example of use:*

    * `Distribution setting regular h (array) <Examples/Geom/distrib1.py>`_:

    .. literalinclude:: ../build/Examples/Geom/distrib1.py

    * `Distribution setting regular h (pyTree) <Examples/Geom/distrib1PT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/distrib1PT.py


---------------------------------------

.. py:function:: Geom.distrib2(a, h1, h2, add=20, forceAdd=False, normalized=True, algo=0)

    Return the distribution enforcing h1 and h2 at curve extremities.

    :param a: 1D mesh
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones] 
    :param h1: size to enforce at first extremity
    :type h1: float
    :param h2: size to enforce at second extremity
    :type h2: float
    :param add: when algo=0, number of points to add starting from regular distribution
    :type add: integer
    :param forceAdd: when algo=0, force add to be exact
    :type forceAdd: boolean
    :param normalized: if True, return a normalized distribution
    :type normalized: boolean
    :param algo: 0: hyperbolic tangent, 1: geometric
    :type algo: integer
    :return: the distribution of the curve as mesh coordinates with h enforced
    :rtype: Identical to a

    *Example of use:*

    * `Distribution setting h1 and h2 (array) <Examples/Geom/distrib2.py>`_:

    .. literalinclude:: ../build/Examples/Geom/distrib2.py

    * `Distribution setting h1 and h2 (pyTree) <Examples/Geom/distrib2PT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/distrib2PT.py


---------------------------------------

.. py:function:: Geom.getDistribution(a)
  
    Return the distribution (curvilinear abscissa) of a curve as a mesh coordinates.
    
    :param a: 1D mesh
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones] 
    :return: the distribution of the curve as mesh coordinates
    :rtype: Identical to a

    *Example of use:*

    * `Distribution of a uniform NACA0012 profile (array) <Examples/Geom/getDistribution.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getDistribution.py

    * `Distribution of a uniform NACA0012 profile (pyTree) <Examples/Geom/getDistributionPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getDistributionPT.py


---------------------------------------------------------------------------

Information about geometries
--------------------------------

For pyTrees, the information is stored as a son node of 'FlowSolution' if it is defined for
all the points of the geometry.


.. py:function:: Geom.getLength(a)
  
    Return the length of a discretized curve or a set of curves.

    :param a: curve or list of curves
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :return: the length of curves
    :rtype: float

    *Example of use:*

    * `Length of a line  (array) <Examples/Geom/getLength.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getLength.py

    * `Length of a line (pyTree) <Examples/Geom/getLengthPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getLengthPT.py

---------------------------------------

.. py:function:: Geom.getDistantIndex(a, ind, l)
  
    Return the point index in a that is distant of l from a point of index ind in a.

    :param a: 1D mesh 
    :type a: array or zone 
    :param ind: index of starting point 
    :type ind: integer
    :param l: distance of the end point to the starting point of index ind
    :type l: float
    :return: the index in a of the end point at distance l of ind
    :rtype: integer

    *Example of use:*

    * `Index of point distant to another one  (array) <Examples/Geom/getDistantIndex.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getDistantIndex.py

    * `Index of point distant to another one (pyTree) <Examples/Geom/getDistantIndexPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getDistantIndexPT.py

---------------------------------------


.. py:function:: Geom.getNearestPointIndex(a, P)
  
    Return the index and the squared distance of the nearest point of P(x,y,z) in a
    If a is a list of meshes, the minimum distance for all meshes in a from P is returned.

    :param a: 1D mesh 
    :type a: [array, list of arrays] or [pyTree,base, list of zones, zone] 
    :param P: coordinates of the point P or point list P 
    :type P: (float,float,float) or [(float,float,float),...,(float,float,float)] 
    :return: the index and squared distance of the nearest point(s) of a to point(s) P
    :rtype: [(integer,float) or list of (integer,float)]

    *Example of use:*

    * `Index of nearest point to P (array) <Examples/Geom/getNearestPointIndex.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getNearestPointIndex.py

    * `Index of nearest point to P (pyTree) <Examples/Geom/getNearestPointIndexPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getNearestPointIndexPT.py

---------------------------------------

.. py:function:: Geom.getCurvatureRadius(a)
  
    Return the curvature radius of a curve a. 

    Exists also as in place version (_getCurvatureRadius) that modifies a and returns None.

    :param a: 1D mesh 
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :return: the curvature radius named 'radius'.
    :rtype: Identical to a

    *Example of use:*

    * `Curvature radius of a curve (array) <Examples/Geom/getCurvatureRadius.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvatureRadius.py

    * `Curvature radius of a curve (pyTree) <Examples/Geom/getCurvatureRadiusPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvatureRadiusPT.py

---------------------------------------

.. py:function:: Geom.getCurvatureAngle(a)
  
    Return the curvature angle of a curve a. 
    
    Exists also as in place version (_getCurvatureAngle) that modifies a and returns None.

    :param a: 1D mesh 
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :return: the curvature angle named 'angle'.
    :rtype: Identical to a

    *Example of use:*

    * `Curvature angle of a curve (array) <Examples/Geom/getCurvatureAngle.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvatureAngle.py

    * `Curvature angle of a curve (pyTree) <Examples/Geom/getCurvatureAnglePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvatureAnglePT.py

---------------------------------------

.. py:function:: Geom.getCurvatureHeight(a)
  
    Return the curvature height of a curve a. 
    
    Exists also as in place version (_getCurvatureHeight) that modifies a and returns None.

    :param a: 1D mesh 
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :return: the curvature height named 'hmax' as an array or as a flow solution at nodes.
    :rtype: Identical to a

    *Example of use:*

    * `Curvature height of a curve (array) <Examples/Geom/getCurvatureHeight.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvatureHeight.py

    * `Curvature height of a curve (pyTree) <Examples/Geom/getCurvatureHeightPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvatureHeightPT.py

---------------------------------------

.. py:function:: Geom.getSharpestAngle(a)
  
    Return the sharpest angle (in degrees) of a curve. Sharpest angle is defined
    at each node of input curve.

    Exists also as in place version (_getSharpestAngle) that modifies a and returns None.

    :param a: 1D mesh 
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones] 
    :return: the sharpest angle named 'alpha' as an array or as a flow solution at nodes.
    :rtype: Identical to a

    *Example of use:*

    * `Sharpest angle of a curve (array) <Examples/Geom/getSharpestAngle.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getSharpestAngle.py

    * `Sharpest angle of a curve (pyTree) <Examples/Geom/getSharpestAnglePT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getSharpestAnglePT.py

---------------------------------------

.. py:function:: Geom.getCurvilinearAbscissa(a)
  
    Return the curvilinear abscissa of a curve a (scalar in range [0.,1.]).
    
    Exists also as in place version (_getCurvilinearAbcissa) that modifies a and returns None.

    :param a: 1D mesh 
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones] 
    :return: the curvilinear abscissa named 's' as an array or as a flow solution at nodes.
    :rtype: Identical to a

    *Example of use:*

    * `Curvilinear abscissa of a curve (array) <Examples/Geom/getCurvilinearAbscissa.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvilinearAbscissa.py

    * `Curvilinear abscissa of a curve (pyTree) <Examples/Geom/getCurvilinearAbscissaPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getCurvilinearAbscissaPT.py

---------------------------------------

.. py:function:: Geom.getTangent(a)
  
    Return the unit tangent vector of all nodes of a 1D array (only structured) as a mesh coordinates.

    :param a: 1D structured mesh
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones] 
    :return: the unit tangent vector of the curve as mesh coordinates
    :rtype: Identical to a

    *Example of use:*

    * `Unit tangent vector of a spline (array) <Examples/Geom/getTangent.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getTangent.py

    * `Unit tangent vector of a spline (pyTree) <Examples/Geom/getTangentPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getTangentPT.py

---------------------------------------

.. py:function:: Geom.getUV(a, normalDeviationWeight=2., texResolution=1920)
  
    Return the UV map of a TRI surface. Return also color and bumpmap atlas as Zone image.
    The number of points of returned zone can be greater than the one of input zone because
    of seams in the map (lines where UV are not unique).

    :param a: single TRI zone
    :type a: [array] or [zone]
    :param normalDeviationWeight: weight for normal deviation to form a chart
    :type normalDeviationWeight: float
    :param texResolution: approx resolution of textures
    :type texResolution: float
    :return: The UV map of a, color atlas, bumpmap atlas
    :rtype: three zones

    *Example of use:*

    * `Get UV map of TRI surface (array) <Examples/Geom/getUV.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getUV.py

    * `Get UV map of TRI surface (pyTree) <Examples/Geom/getUVPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getUVPT.py

---------------------------------------

.. py:function:: Geom.getUVFromIJ(a)
  
    Return the UV map of a STRUCT surface. UV is based on IJ of surface.

    :param a: input
    :type a: [array, arrays] or [zone, zones, base, tree]
    :rtype: identical to input

    *Example of use:*

    * `Get UV map from IJ (array) <Examples/Geom/getUVFromIJ.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getUVFromIJ.py

    * `Get UV map from IJ (pyTree) <Examples/Geom/getUVFromIJPT.py>`_:

    .. literalinclude:: ../build/Examples/Geom/getUVFromIJPT.py

---------------------------------------

.. toctree::
   :maxdepth: 2   


Index
#######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


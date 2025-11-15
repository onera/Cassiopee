# - addSpline (array) -
import OCC
import numpy

# control points
points = numpy.array([[0,0,0],[1,1,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0]], dtype=numpy.float64, order='F')

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")

# polyline
OCC.occ.addSpline(hook, points, 0, 1)

# from control points with chord length param
OCC.occ.addSpline(hook, points, 0, 2)

# from control points with uniform param
OCC.occ.addSpline(hook, points, 2, 3)

# interpolation from discrete 1D mesh
import Geom as D
a = D.naca(12., 50)
OCC.occ.addSpline(hook, a[1], 1, 3)

OCC.occ.writeCAD(hook, "out.step", "fmt_step")

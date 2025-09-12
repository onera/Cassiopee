# - addSpline (array) -
import OCC
import numpy

# from a list of points
hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
points = numpy.array([[0,0,0],[1,1,0],[2,0,0],[3,0,0]], dtype=numpy.float64, order='F')
OCC.occ.addSpline(hook, points, 1)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

# from a discrete 1D mesh
import Geom as D
a = D.naca(12., 50)
OCC.occ.addSpline(hook, a[1], 1)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

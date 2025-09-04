# - addSpline (array) -
import OCC
import numpy

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
knots = numpy.array([[0,0,0],[1,1,0],[2,0,0],[3,0,0]], dtype=numpy.float64, order='F')
OCC.occ.addSpline(hook, knots, 1)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

# - addSpline (array) -
import OCC
import numpy

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
knots = numpy.array([[0,0,0],[1,1,0],[2,0,0]], dtype=numpy.float64)
OCC.occ.addSpline(hook, knots)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

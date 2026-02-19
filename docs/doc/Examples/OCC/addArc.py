# - addArc (array) -
import OCC
import numpy
h = numpy.sqrt(2.)/2.
hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addArc(hook, (0,0,0), (h,1-h,0), (1,1,0))
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

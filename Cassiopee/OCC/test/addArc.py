# - addArc (array) -
import OCC
import numpy
h = numpy.sqrt(2.)/2.
hook = OCC.createEmptyCAD()
OCC._addArc(hook, (0,0,0), (h,1-h,0), (1,1,0))
OCC.writeCAD(hook, "out.step")

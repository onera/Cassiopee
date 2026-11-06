# - revolve (array) -
import OCC
import numpy
h = numpy.sqrt(2.)/2.

hook = OCC.createEmptyCAD()
OCC._addArc(hook, (-1,0,0), (-h,h,0), (0,1,0))
OCC._addLine(hook, (0,1,0), (1,1,0))

OCC._revolve(hook, [1,2], (0,0,0), (1,0,0), 360.)

OCC.writeCAD(hook, "out.step")

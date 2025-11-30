# - revolve (array) -
import OCC
import numpy
h = numpy.sqrt(2.)/2.

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addArc(hook, (-1,0,0), (-h,h,0), (0,1,0))
OCC.occ.addLine(hook, (0,1,0), (1,1,0))

OCC.occ.revolve(hook, [1,2], (0,0,0), (1,0,0), 360.)

OCC.occ.writeCAD(hook, "out.step", "fmt_step")

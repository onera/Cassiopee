# - splitEdge (array) -
import OCC
import Geom as D
import numpy

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")

a = D.naca(12., 50)
OCC.occ.addSpline(hook, a[1], 1, 3)

OCC.occ.splitEdge(hook, 1, -999., 0.0, 0.0, 0.0)

OCC.occ.writeCAD(hook, "out.step", "fmt_step")

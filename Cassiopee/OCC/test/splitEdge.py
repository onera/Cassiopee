# - splitEdge (array) -
import OCC
import Geom as D

hook = OCC.createEmptyCAD()

a = D.naca(12., 50)
OCC._addSpline(hook, a[1], 1, 3)

OCC._splitEdge(hook, 1, -999., (0.0,0.0,0.0))

OCC.writeCAD(hook, "out.step")

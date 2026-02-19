# - addCircle (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)
OCC.writeCAD(hook, "out.step")

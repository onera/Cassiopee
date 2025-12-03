# - fillHole (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)
OCC._fillHole(hook, [1], [], 0)
OCC.writeCAD(hook, "out.step")

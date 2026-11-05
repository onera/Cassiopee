# - reverse (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addBox(hook, (0,0,0), 1., 1., 1.)
OCC._reverse(hook)
OCC.writeCAD(hook, "out.step")

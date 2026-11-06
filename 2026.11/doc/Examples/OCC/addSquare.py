# - addSquare (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addSquare(hook, (0,0,0), 1., 1.)
OCC.writeCAD(hook, "out.step")

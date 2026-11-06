# - addSquare2 (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addSquare2(hook, (0,0,0), (1,0,0), (1,1,0), (0,1,0))
OCC.writeCAD(hook, "out.step")

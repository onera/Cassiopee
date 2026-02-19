# - addLine (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (1,0,0))
OCC.writeCAD(hook, "out.step")

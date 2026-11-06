# - addCylinder (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addCylinder(hook, (0,0,0), (0,0,1), 1., 3.)
OCC.writeCAD(hook, "out.step")

# - bottle (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC.occ.bottle(hook, 1., 1., 0.5)
OCC.writeCAD(hook, "out.step")

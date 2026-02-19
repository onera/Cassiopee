# - addCylinder (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addCylinder(hook, (0,0,0), (0,0,1), 1, 3)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

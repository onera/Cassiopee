# - bottle (array) -
import OCC

hook = OCC.occ.createEmptyCAD("bottle.stp", "fmt_step")
OCC.occ.bottle(hook, 1., 1., 0.5)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

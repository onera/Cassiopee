# - addSphere (array) -
import OCC

hook = OCC.occ.createEmptyCAD("bottle.stp", "fmt_step")
OCC.occ.addSphere(hook, 0., 0., 0., 1.)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

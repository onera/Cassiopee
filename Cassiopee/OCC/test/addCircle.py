# - addCircle (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addCircle(hook, (0,0,0), (0,0,1), 1., False)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

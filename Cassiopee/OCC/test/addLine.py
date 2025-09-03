# - addLine (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addLine(hook, (0,0,0), (1,0,0))
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

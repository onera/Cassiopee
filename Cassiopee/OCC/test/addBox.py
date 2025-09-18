# - addBox (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addBox(hook, (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1))
OCC.occ.writeCAD(hook, "out.step", "fmt_step")

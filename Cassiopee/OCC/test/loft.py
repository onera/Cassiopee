# - loft (array) -
import OCC

hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addCircle(hook, (0,0,0), (0,0,1), 1., False)
OCC.occ.addCircle(hook, (0,0,1), (0,0,1), 1., False)
OCC.occ.addCircle(hook, (0,0,2), (0,0.3,0.5), 1., False)

OCC.occ.loft(hook, [1,2,3], [])

OCC.occ.writeCAD(hook, "out.step", "fmt_step")

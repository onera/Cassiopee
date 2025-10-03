# - loft (array) -
import OCC

N = 10
out = []
for i in range(N):
    hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
    OCC.occ.addCircle(hook, (0,0,0), (0,0,1), 1., False)
    OCC._rotate(hook, (20,0,0), (0,1,0), i/N*90.)
    out.append(hook)

hook = OCC.occ.mergeCAD(out)
OCC.occ.loft(hook, [i for i in range(1,N)], [])

OCC.occ.writeCAD(hook, "out.step", "fmt_step")

# - intersectEdgeFace (array) -
import OCC
import numpy
from math import sin, cos

# sketches
N = 5
ucurves = []
for i in range(N):
    hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
    OCC.occ.addCircle(hook, (20,0,0), (0,0,1), 1., False)
    OCC._rotate(hook, (0,0,0), (0,1,0), i/(N-1)*90.)
    ucurves.append(hook)

# y = 0 plane
hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
OCC.occ.addSquare(hook, (-2,0,-22), (22,0,-22), (22,0,2), (-2,0,2), True)

hook = OCC.occ.mergeCAD(ucurves+[hook])
for i in range(N):
    ret = OCC.occ.intersectEdgeFace(hook, [i+1], [1])
    print(ret)

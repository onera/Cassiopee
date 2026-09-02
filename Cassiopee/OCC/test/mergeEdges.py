# - mergeEdges (array) -
import OCC

# Figure with edge 1 and 2 to merge
hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (0.5,0,0))
OCC._addLine(hook, (0.5,0,0), (1.,0,0))
OCC._addLine(hook, (1,0,0), (1,1,0))
OCC._addLine(hook, (1,1,0), (0,1,0))
OCC._addLine(hook, (0,1,0), (0,0,0))
OCC._fillHole(hook, [1,2,3,4,5], [], 0)

OCC._mergeEdges(hook, [1,2])

OCC.writeCAD(hook, "out.step")

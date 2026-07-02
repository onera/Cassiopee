# - removeEdges (array) -
import OCC

hook = OCC.createEmptyCAD()

OCC._addLine(hook, (0,0,0), (1,0,0))
OCC._addLine(hook, (1,0,0), (1,1,0))
OCC._addLine(hook, (1,1,0), (0,1,0))

OCC._removeEdges(hook, [2])

OCC.writeCAD(hook, 'out.step')

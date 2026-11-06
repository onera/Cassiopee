# - booleanEdges (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)
OCC._addCircle(hook, (1,0,0), (0,0,1), 0.7)
ret = OCC._booleanEdges(hook, [1], [2], op=0)
OCC.writeCAD(hook, 'out.step')

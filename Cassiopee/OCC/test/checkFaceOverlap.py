# - checkFaceOverlap (array) -
import OCC

hook = OCC.createEmptyCAD()

OCC._addSquare(hook, (0,0,0), 1., 1., makeFace=True)
OCC._addSquare(hook, (0,0,0), 1., 1., makeFace=True)

OCC.checkFaceOverlap(hook, 1, 2)

OCC.writeCAD(hook, 'out.step')

# - meshAllOCC (array) -
import Converter as C
import OCC

hook = OCC.createEmptyCAD()
OCC._addBox(hook, (0,0,0), 1., 1., 1.)
m = OCC.meshAllOCC(hook, 0.1)
C.convertArrays2File(m[1], 'out.plt')

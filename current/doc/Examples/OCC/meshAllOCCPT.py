# - meshAllOCC (pyTree) -
import Converter.PyTree as C
import OCC.PyTree as OCC

hook = OCC.createEmptyCAD()
OCC._addBox(hook, (0,0,0), 1., 1., 1.)
m = OCC.meshAllOCC(hook, 0.1)
C.convertPyTree2File(m, 'out.cgns')

# - getFaceNameInOCAF (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C
#import OCC.occ as occ

hook = OCC.readCAD("cube.step", "fmt_step")
#occ.printOCAF(hook)

#ret = occ.getFaceNameInOCAF(hook)
#print(ret)
t = OCC.meshAll(hook, hmax=2.)
OCC._addOCAFCompoundNames(hook, t)
C.convertPyTree2File(t, 'out.cgns')
# - getFaceNameInOCAF (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C

hook = OCC.readCAD("cube.step", "fmt_step")
#OCC.printOCAF(hook)

ret = OCC.getFaceNameInOCAF(hook)
print(ret)

t = OCC.meshAll(hook, hmax=2.)
OCC._addOCAFCompoundNames(hook, t)
C.convertPyTree2File(t, 'out.cgns')

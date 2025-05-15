# - meshAll (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C

hook = OCC.readCAD("cube.step", "fmt_step")
# constant h
m = OCC.meshAll(hook, hmin=1., hmax=1.)
# varying h between hmin and hmax with deflection hausd
m = OCC.meshAll(hook, hmin=1., hmax=5., hausd=1.)
C.convertPyTree2File(m, 'out.cgns')


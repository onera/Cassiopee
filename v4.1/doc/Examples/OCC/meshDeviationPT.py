# - meshDeviation (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C

hook = OCC.readCAD("cube.step", "fmt_step")
m = OCC.meshAll(hook, hmin=1., hmax=1.)

OCC._meshDeviation(hook, m)

C.convertPyTree2File(m, 'out.cgns')

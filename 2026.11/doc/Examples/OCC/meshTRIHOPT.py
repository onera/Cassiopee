# - meshTRIHO (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C

m = OCC.meshTRIHO("cube.step", "fmt_step")
C.convertPyTree2File(m, 'out.cgns')

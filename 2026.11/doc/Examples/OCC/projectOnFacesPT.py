# - projectOnFaces (PyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C
import Generator.PyTree as G

hook = OCC.readCAD("cube.step", "fmt_step")

# project on all CAD faces
a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC._projectOnFaces(hook, a, None)
C.convertPyTree2File(a, 'out.cgns')

# project on certain CAD faces
a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC._projectOnFaces(hook, a, [1,2,3,4,5,6])
C.convertPyTree2File(a, 'out.cgns')

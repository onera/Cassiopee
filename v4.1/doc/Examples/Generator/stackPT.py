# - stack (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T

a = G.cylinder((0,0,0), 1, 1.3, 360, 0, 1., (50,10,1))
b = T.rotate(a, (0,0,0), (1,0,0), 5.)
b = T.translate(b, (0,0,0.5))
c = G.stack(a, b)
C.convertPyTree2File(c, 'out.cgns')

# - splitTRI (PyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T

a = D.circle( (0,0,0), 1, N=20 )
a = C.convertArray2Tetra(a)
a = G.close(a)
b = G.T3mesher2D(a)
c = [[9, 25, 27, 30, 29, 28, 34, 38, 0], [29, 23, 19, 20, 24, 29]]
D = T.splitTRI(b, c)
C.convertPyTree2File(D, 'out.cgns')

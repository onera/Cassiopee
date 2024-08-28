# - projectOrtho (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = D.sphere((0,0,0), 1., 20)
b = G.cart((1.1,-0.1,-0.1),(0.1,0.1,0.1), (1,5,5))
c = T.projectOrtho(b, a); c[0] = 'projection'
C.convertPyTree2File([a,b,c], 'out.cgns')

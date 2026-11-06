import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Post.PyTree as P

a = D.naca(12., N=101)
b = G.cart((-15,-15,0), (1,1,1), (31,31,1))
b = P.exteriorFacesStructured(b)
C.convertPyTree2File([a]+b, 'try1DGeometry.cgns')

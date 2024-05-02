# - exteriorFacesStructured (pyTree)-
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (4,4,6))
zones = P.exteriorFacesStructured(a)
C.convertPyTree2File(zones, 'out.cgns')

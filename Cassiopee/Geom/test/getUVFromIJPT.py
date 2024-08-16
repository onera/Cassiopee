# - getUVFromIJ (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))

D._getUVFromIJ(a)

C.convertPyTree2File(a, 'out.cgns')

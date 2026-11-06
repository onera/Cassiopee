# - splitNParts (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (81,81,81))
b = G.cart((80,0,0), (1,1,1), (41,81,41))
t = C.newPyTree(['Base',a,b])
t = T.splitNParts(t, 10, multigrid=0, dirs=[1,2,3])
C.convertPyTree2File(t, 'out.cgns')

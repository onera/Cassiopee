# - splitFullMatch (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Connector.PyTree as X

z0 = G.cart((0.,0.,0.), (1,1,1), (10,10,1))
z1 = G.cart((0.,9.,0.),(1,1,1),(5,10,1))
z2 = G.cart((4.,9.,0.),(1,1,1),(6,10,1))
t = C.newPyTree(['Base',z0,z1,z2])
t = X.connectMatch(t, dim=2)

T._splitFullMatch(t)

C.convertPyTree2File(t, 'out.cgns')

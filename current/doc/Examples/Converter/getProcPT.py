# - getProc (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')

zones = Internal.getNodesFromType(t, 'Zone_t')
for z in zones: print(z[0]+' -> '+str(Cmpi.getProc(z)))
#>> cart -> 0
#>> cart.0 -> 1

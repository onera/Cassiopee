# - setProc (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
t = Cmpi.setProc(t, 1)

zones = Internal.getZones(t)
for z in zones: print(z[0]+' -> '+str(Cmpi.getProc(z)))
#>> cart -> 1
#>> cart.0 -> 1

# - writeZones (pyTree) -
import Converter.PyTree as C
import Converter.Distributed as Distributed
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base'])
C.convertPyTree2File(t, 'out.adf')
t[2][1][2] += [a,b]
(t, dic) = Distributor2.distribute(t, NProc=2, algorithm='fast')
Distributed.writeZones(t, 'out.adf', proc=0)

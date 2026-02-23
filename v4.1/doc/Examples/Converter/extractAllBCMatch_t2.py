# - extractAllBCMatch (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

# NGON
a = G.cartNGon((1,1,1), (1.,1.,1.), (4,10,3)); a[0]='cart1'
b = G.cartNGon((4,2,0), (1.,1.,1.), (5, 8,5)); b[0]='cart2'

t = C.newPyTree(['Base',a,b])

t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = C.initVars(t, '{centers:G}=2.53')
t = C.initVars(t, '{centers:H}={centers:CoordinateY}')
t = C.initVars(t, '{centers:M}={centers:CoordinateX}')
t = X.connectMatch(t,dim=3)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall')

dico = C.extractAllBCMatch(t,['centers:H','centers:M'])
test.testO(dico, 1)

dico = C.extractAllBCMatch(t)
test.testO(dico, 2)

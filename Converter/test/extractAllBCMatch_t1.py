# - extractAllBCMatch (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as CP
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

a = G.cart((1,1,1), (1.,1.,1.), (4,10,3)); a[0]='cart1'
b = G.cart((4,2,0), (1.,1.,1.), (5, 8,5)); b[0]='cart2'

a = T.reorder(a,(-3,1,-2))

t = CP.newPyTree(['Base',a,b])

t = CP.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = CP.initVars(t, '{centers:G}=2.53')
t = CP.initVars(t, '{centers:H}={centers:CoordinateY}')
t = CP.initVars(t, '{centers:M}={centers:CoordinateX}')
t = X.connectMatch(t,dim=3)
t = CP.fillEmptyBCWith(t,"wall",'BCWall')

dico = CP.extractAllBCMatch(t,['centers:G','centers:H','centers:M']) 

# print dico

test.testO(dico, 1)


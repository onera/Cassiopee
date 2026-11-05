# - extractAllBCMatch (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import KCore.test as test

varL = ['centers:H','centers:M']

a = G.cartNGon((1,1,1), (1.,1.,1.), (4,10,3)); a[0]='cart1'
b = G.cartNGon((4,2,0), (1.,1.,1.), (5, 8,5)); b[0]='cart2'

t = C.newPyTree(['Base',a,b])

t = C.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = C.initVars(t, '{centers:G}=2.3')
t = C.initVars(t, '{centers:H}={centers:CoordinateY}')
t = C.initVars(t, '{centers:M}={centers:CoordinateX}')
t = X.connectMatch(t, dim=3)
t = C.fillEmptyBCWith(t, "wall", 'BCWall')

dico = C.extractAllBCMatch(t, varL)
it   = 0

for z in Internal.getZones(t):
    it  += 1
    indR, fld = C.computeBCMatchField(z,dico,varL)

    test.testO(indR, it)
    test.testO(fld,  it+2)

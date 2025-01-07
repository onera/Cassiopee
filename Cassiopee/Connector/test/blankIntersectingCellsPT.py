# - blankIntersectingCells (pyTree)
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Connector.PyTree as X

a1 = G.cart((0.,0.,0.),(1.,1.,1.),(11,11,11))
a2 = T.rotate(a1, (0.,0.,0.), (0.,0.,1.),10.)
a2 = T.translate(a2, (7.,5.,5.)); a1[0] = 'cart1'; a2[0] = 'cart2'
t = C.newPyTree(['Base',a1,a2])
C._initVars(t,'centers:cellN',1.)
t2 = X.blankIntersectingCells(t, tol=1.e-10)
C.convertPyTree2File(t2, "out.cgns")

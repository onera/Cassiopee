# - blankIntersectingCells (array)
import Converter as C
import Generator as G
import Transform as T
import Connector as X
a1 = G.cart((0.,0.,0.),(1.,1.,1.),(11,11,11))
a2 = T.rotate(a1, (0.,0.,0.), (0.,0.,1.),10.)
a2 = T.translate(a2, (7.,5.,5.))
A = [a1,a2]
Ac = C.node2Center(A); Ac = C.initVars(Ac,'cellN',1.);
Ac = X.blankIntersectingCells(A, Ac, tol=1.e-10)
C.convertArrays2File(Ac,"out.plt")

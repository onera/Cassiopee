# - getMaxLength(pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as T

a = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
a = G.enforcePlusX(a,1e-6,(5,50))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',1.)
a = G.getMaxLength(a)
T.testT(a,1)
#
# TEST NON STRUCTURE ELT BASIQUE
a = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
a = G.enforcePlusX(a,1e-6,(5,50))
a = C.convertArray2Tetra(a)
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',1.)
a = G.getMaxLength(a)
T.testT(a,2)
#
# TEST NON STRUCTURE NGON
a = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
a = G.enforcePlusX(a,1e-6,(5,50))
a = C.convertArray2NGon(a)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',1.)
a = G.getMaxLength(a)
T.testT(a,3)

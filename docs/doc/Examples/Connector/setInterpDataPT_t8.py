# - setInterpData (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

# Donor
zD = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
C._initVars(zD,'centers:cellN',1.)

# Receptor
zR = G.cart((0.25,0.25,0.25), (0.1,0.1,0.1), (5,5,5))
C._initVars(zR,'{centers:cellN}=2.*({centers:CoordinateX}<1.)')
zR=X.setInterpData(zR, zD, order=2, penalty=1, nature=0, method='conservative', hook=None, dim=3, loc='centers',itype='chimera')
test.testT(zR,1)

zD = X.setInterpData(zR, zD, order=2, penalty=1, nature=0, method='conservative', hook=None, dim=3, loc='centers',itype='chimera',storage='inverse')
test.testT(zD,2)

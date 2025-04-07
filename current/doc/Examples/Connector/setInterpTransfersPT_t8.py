import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

# 1. Direct storage
# Donor
zD = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
C._initVars(zD,'centers:cellN',1.)

# Receptor
zR = G.cart((0.25,0.25,0.25), (0.1,0.1,0.1), (5,5,5))
C._initVars(zR,'{centers:cellN}=2.*({centers:CoordinateX}<1.)')
zR=X.setInterpData(zR, zD, method='conservative',itype='chimera')
C._initVars(zR,"centers:F", 0.)
C._initVars(zD,'{centers:F}={centers:CoordinateX}*{centers:CoordinateY}')
zD = C.node2Center(zD)
zR2 = X.setInterpTransfers(zR,zD,variables=["F"])
test.testT(zR2,1)
# in place
X._setInterpTransfers(zR,zD,variables=["F"])
test.testT(zR,2)


# Inverse storage
# Donor
zD = G.cart((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
C._initVars(zD,'centers:cellN',1.)

# Receptor
zR = G.cart((0.25,0.25,0.25), (0.1,0.1,0.1), (5,5,5))
C._initVars(zR,'{centers:cellN}=2.*({centers:CoordinateX}<1.)')
zD=X.setInterpData(zR, zD, method='conservative',itype='chimera',storage='inverse')
C._initVars(zR,"centers:F", 0.)
C._initVars(zD,'{centers:F}={centers:CoordinateX}*{centers:CoordinateY}')
zD = C.node2Center(zD)
zR2 = X.setInterpTransfers(zR,zD,variables=["F"])
test.testT(zR2,3)
# in place
X._setInterpTransfers(zR,zD,variables=["F"])
test.testT(zR,4)

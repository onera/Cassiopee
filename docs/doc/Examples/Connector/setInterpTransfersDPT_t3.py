# - setIBCTransfersD (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Post.PyTree as P
import numpy as N
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((-1,-1,-1),(0.04,0.04,1),(51,51,3))
s = G.cylinder((0,0,-1), 0, 0.4, 360, 0, 4, (15,15,5))
s = C.convertArray2Tetra(s); s = T.join(s); s = P.exteriorFaces(s)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
# Blanking
bodies = [[s]]
BM = N.array([[1]], Internal.E_NpyInt)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
X._setHoleInterpolatedPoints(t,depth=-2)

# Dist2Walls
DTW._distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
C._rmVars(t,['TurbulentDistance'])
C._initVars(t,"centers:Density",1.)
C._initVars(t,"centers:VelocityX",0.2)
C._initVars(t,"centers:VelocityY",0.)
C._initVars(t,"centers:VelocityZ",0.)
C._initVars(t,"centers:Temperature",1.)
tc = C.node2Center(t)
X._setIBCData(t, tc, loc='centers', storage='inverse')
vars = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']

tc2 = Internal.copyRef(tc)
zones = Internal.getNodesFromType2(t, 'Zone_t')
X.miseAPlatDonorTree__(zones, tc2, graph=None)
info = X.setInterpTransfersD(tc2, bcType=0, varType=2,variablesIBC=vars,compact=0, compactD=1)
test.testO(info)
test.testA([info[0][1]],2)
#
# variable turbulente SA
#
vars+=['TurbulentSANuTilde']
C._initVars(tc,'TurbulentSANuTilde',15.)
C._initVars(t,"centers:TurbulentSANuTilde",15.)
zones = Internal.getNodesFromType2(t, 'Zone_t')
X.miseAPlatDonorTree__(zones, tc, graph=None)
info = X.setInterpTransfersD(tc, bcType=0,varType=21,variablesIBC=vars, compact=0, compactD=1)
test.testO(info,3)
test.testA([info[0][1]],4)
C._initVars(tc,'TurbulentSANuTilde',15.)

info = X.setInterpTransfersD(tc, bcType=1,varType=21,variablesIBC=vars, compact=0, compactD=1)
test.testO(info,5)
test.testA([info[0][1]],6)

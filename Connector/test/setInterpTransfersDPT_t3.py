# - setIBCTransfersD (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Post.PyTree as P
import numpy as N
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import KCore.test as test

a = G.cart((-1,-1,-1),(0.04,0.04,1),(51,51,3))
s = G.cylinder((0,0,-1), 0, 0.4, 360, 0, 4, (15,15,5)) 
s = C.convertArray2Tetra(s); s = T.join(s); s = P.exteriorFaces(s)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
# Blanking
bodies = [[s]]
BM = N.array([[1]],N.int32)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
t = X.setHoleInterpolatedPoints(t,depth=-2)

# Dist2Walls
t = DTW.distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
t = C.rmVars(t,['TurbulentDistance'])
t = I.initConst(t,MInf=0.2,loc='centers')
tc = C.node2Center(t)
tc = X.setIBCData(t, tc, loc='centers', storage='inverse')
info = X.setInterpTransfersD(tc)
test.testO(info)
test.testA([info[0][1]],2)
#
# variable turbulente SA
#
tc = C.initVars(tc,'TurbulentSANuTilde',15.)
vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity', 'TurbulentSANuTilde']
info = X.setInterpTransfersD(tc, bcType=0,varType=11,variablesIBC=vars)
test.testO(info,3)
test.testA([info[0][1]],4)

info = X.setInterpTransfersD(tc, bcType=1,varType=11,variablesIBC=vars)
test.testO(info,5)
test.testA([info[0][1]],6)

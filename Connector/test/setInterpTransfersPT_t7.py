# - setInterpTransfers IBC + extraction (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Post.PyTree as P
import numpy as N
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((-1,-1,-1),(0.04,0.04,1),(51,51,3))
s = G.cylinder((0,0,-1), 0, 0.4, 360, 0, 4, (30,30,5)) 
s = C.convertArray2Tetra(s); s = T.join(s); s = P.exteriorFaces(s)
t = C.newPyTree(['Base',a])
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
t = I.initConst(t,MInf=0.2,loc='centers')
tc = C.node2Center(t)
no = 1
for bcType in range(4):
    tp = X.setIBCData(t, tc, loc='centers', storage='direct', bcType=bcType)
    for varType in range(1,4):
        t2 = X.setInterpTransfers(tp, tc, bcType=bcType, varType=varType)
        test.testT(t2,no)
        no+=1

# Stockage inverse
Internal._rmNodesByName(t,"IBCD_*")
for bcType in range(4):
    tc2 = X.setIBCData(t, tc, loc='centers', storage='inverse', bcType=bcType)
    for varType in range(1,4):
        X._setInterpTransfers(t, tc2, bcType=bcType, varType=varType)
        test.testT(tc, no)
        no += 1

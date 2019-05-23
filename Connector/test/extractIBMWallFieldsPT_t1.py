# - extractionIBM a la paroi (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Converter.Internal as Internal
import Connector.ToolboxIBM as IBM
import KCore.test as test
import numpy , sys

N = 41
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,N))
xm = 0.5*N/(N-1)
s = D.sphere((xm,xm,xm),0.1,N=20)
s = C.convertArray2Tetra(s); s = G.close(s)
t = C.newPyTree(['Base', a])

# Blanking
bodies = [[s]]
BM = numpy.array([[1]],numpy.int32)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
t = X.setHoleInterpolatedPoints(t,depth=-2)
# Dist2Walls
DTW._distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
I._initConst(t,MInf=0.2,loc='centers')
tc = C.node2Center(t)

tb = C.newPyTree(['Base', s])
C._addState(tb, 'EquationDimension',3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')

tp = X.setIBCData(t, tc, loc='centers', storage='direct', bcType=0)
t2 = X.setInterpTransfers(tp, tc, bcType=0, varType=1)
z = IBM.extractIBMWallFields(t2, tb=tb)
test.testT(z,1)

#
tp = X.setIBCData(t, tc, loc='centers', storage='direct', bcType=3)
t2 = X.setInterpTransfers(tp, tc, bcType=3, varType=1)
z = IBM.extractIBMWallFields(t2, tb=tb)
test.testT(z,2)

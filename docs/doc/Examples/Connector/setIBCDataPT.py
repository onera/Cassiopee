# - setIBCData (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Post.PyTree as P
import numpy as N
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T

a = G.cart((-1,-1,-1),(0.01,0.01,1),(201,201,3))
s = G.cylinder((0,0,-1), 0, 0.4, 360, 0, 4, (30,30,5))
s = C.convertArray2Tetra(s); s = T.join(s); s = P.exteriorFaces(s)
t = C.newPyTree(['Base',a])
# Blanking
bodies = [[s]]
BM = N.array([[1]],N.int32)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
t = X.setHoleInterpolatedPoints(t, depth=-2)
# Dist2Walls
DTW._distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
t = X.setIBCData(t, t, loc='centers', storage='direct',hi=0.03)
C.convertPyTree2File(t, "out.cgns")

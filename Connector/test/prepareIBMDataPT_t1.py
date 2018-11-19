# - prepareIBMData (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.ToolboxIBM as IBM
import Post.PyTree as P
import Geom.PyTree as D
import Dist2Walls.PyTree as DTW
import KCore.test as test
N = 21
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,N))
body = D.sphere((0.5,0,0),0.1,N=20)
t = C.newPyTree(['Base',a])
tb = C.newPyTree(['Base',body])
C._addState(tb, 'EquationDimension',3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
DTW._distance2Walls(t,bodies=tb,loc='centers',type='ortho')
t,tc=IBM.prepareIBMData(t,tb, DEPTH=2,frontType=0)
res = IBM.extractIBMInfo(tc)
test.testT(tc,1)

# front 1
t = C.newPyTree(['Base',a])
tb = C.newPyTree(['Base',body])
C._addState(tb, 'EquationDimension',3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
DTW._distance2Walls(t,bodies=tb,loc='centers',type='ortho')
t,tc=IBM.prepareIBMData(t,tb, DEPTH=2,frontType=1,shiftIBM=False)
test.testT(tc,12)

# front 1 decale
t = C.newPyTree(['Base',a])
tb = C.newPyTree(['Base',body])
C._addState(tb, 'EquationDimension',3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
DTW._distance2Walls(t,bodies=tb,loc='centers',type='ortho')
t,tc=IBM.prepareIBMData(t,tb, DEPTH=2,frontType=1,shiftIBM=True)
test.testT(tc,13)

# CAS 2D
N = 21
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,2))
a[0]='cart.4'
body = D.circle((0.5,0,0),0.1)
t = C.newPyTree(['Base',a])
tb = C.newPyTree(['Base',body])
C._addState(tb, 'EquationDimension',2)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
DTW._distance2Walls(t,bodies=tb,loc='centers',type='ortho',dim=2)
DEPTH = 2; dimPb = 2; model = 'NSTurbulent'
t,tc=IBM.prepareIBMData(t,tb,DEPTH=2,frontType=0)
test.testT(tc,2)
t = C.newPyTree(['Base',a])
DTW._distance2Walls(t,bodies=tb,loc='centers',type='ortho',dim=2)
DEPTH = 2; dimPb = 2; model = 'NSTurbulent'
t,tc=IBM.prepareIBMData(t,tb,DEPTH=2,frontType=1)
test.testT(tc,22)

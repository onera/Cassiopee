# - prepareIBMData_legacy (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.ToolboxIBM as IBM
import Post.PyTree as P
import Geom.PyTree as D
import Dist2Walls.PyTree as DTW
import KCore.test as test
N = 51
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,N))
body = D.sphere((0.5,0,0),0.1,N=50)
t = C.newPyTree(['Base', a])
tb = C.newPyTree(['Base', body])
C._addState(tb, 'EquationDimension', 3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
t,tc=IBM.prepareIBMData_legacy(t, tb, DEPTH=2, frontType=0)
test.testT(tc,1)
# CAS 2D
N = 21
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,2))
body = D.circle((0.5,0,0),0.1)
t = C.newPyTree(['Base', a])
tb = C.newPyTree(['Base', body])
C._addState(tb, 'EquationDimension', 2)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
DEPTH = 2; dimPb = 2; model = 'NSTurbulent'
t,tc=IBM.prepareIBMData_legacy(t, tb, DEPTH=2, frontType=0)
test.testT(tc,2)

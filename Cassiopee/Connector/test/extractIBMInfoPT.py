# - prepareIBMData_legacy (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.ToolboxIBM as IBM
import Post.PyTree as P
import Geom.PyTree as D
import Dist2Walls.PyTree as DTW

N = 21
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,N))
body = D.sphere((0.5,0,0),0.1,N=20)
t = C.newPyTree(['Base',a])
tb = C.newPyTree(['Base',body])
tb = C.addState(tb, 'EquationDimension',3)
tb = C.addState(tb, 'GoverningEquations', 'NSTurbulent')
t = DTW.distance2Walls(t,bodies=tb,loc='centers',type='ortho')
t = P.computeGrad(t,'centers:TurbulentDistance')
t,tc=IBM.prepareIBMData_legacy(t,tb, DEPTH=2)
res = IBM.extractIBMInfo(tc)
C.convertPyTree2File(res,"res.cgns")

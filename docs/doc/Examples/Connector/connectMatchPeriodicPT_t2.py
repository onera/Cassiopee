# - connectMatchPeriodic (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import KCore.test as test

# Un seul bloc
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
a = C.convertArray2NGon(a); a =G.close(a)
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
t = C.newPyTree(['Base',a])
C._addState(t[2][1], 'EquationDimension', 3)
t = X.connectMatchPeriodic(t,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.],unitAngle=None)
test.testT(t,1)

# sur une base
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
a = C.convertArray2NGon(a); a=G.close(a)
b = C.convertArray2NGon(b); b=G.close(b)
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
C._initVars(t,'F',1.); C._initVars(t,'centers:G',2.)
C._addState(t[2][1], 'EquationDimension', 3)
t[2][1] = X.connectMatchPeriodic(t[2][1],rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.])
t[2][1] = X.connectMatchPeriodic(t[2][1],translation=[0,0,5])
test.testT(t,2)

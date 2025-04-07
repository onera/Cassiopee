# - connectMatch (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import KCore.test as test
import math

# Un seul bloc
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']); t[2][1][2] += [a]
C._addState(t[2][1], 'EquationDimension', 3)
t = X.connectMatchPeriodic(t,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.])
test.testT(t,1)
# 2 blocs coincidents
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 30., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']) ; t[2][1][2] += [a]
C._addState(t[2][1], 'EquationDimension', 3)
t=X.connectMatchPeriodic(t,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.])
test.testT(t,11)

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
C._addBC2Zone(b,'wall','BCWall','jmin')
C._addBC2Zone(b,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
C._addState(t[2][1], 'EquationDimension', 3)
t = X.connectMatchPeriodic(t,translation=[0,0,5])
test.testT(t,2)
# sur une liste de zones
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
C._addBC2Zone(b,'wall','BCWall','jmin')
C._addBC2Zone(b,'overlap','BCOverlap','jmax')
[a,b] = X.connectMatchPeriodic([a,b], rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,90.])
[a,b] = X.connectMatchPeriodic([a,b], translation=[0,0,5])
C._addState(t[2][1], 'EquationDimension',3)
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
test.testT(t,4)
#
# sur une base
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
C._addBC2Zone(b,'wall','BCWall','jmin')
C._addBC2Zone(b,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
C._addState(t[2][1], 'EquationDimension', 3)
t[2][1] = X.connectMatchPeriodic(t[2][1],rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.])
t[2][1] = X.connectMatchPeriodic(t[2][1],translation=[0,0,5])
test.testT(t,5)
#
# sur une base en radian
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
C._addBC2Zone(b,'wall','BCWall','jmin')
C._addBC2Zone(b,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
C._addState(t[2][1], 'EquationDimension', 3)
t[2][1] = X.connectMatchPeriodic(t[2][1],rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,math.pi/2],unitAngle='Radian')
t[2][1] = X.connectMatchPeriodic(t[2][1],translation=[0,0,5])
test.testT(t,6)
#
# sur une base en degres
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
C._addBC2Zone(b,'wall','BCWall','jmin')
C._addBC2Zone(b,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
C._addState(t[2][1], 'EquationDimension', 3)
t[2][1] = X.connectMatchPeriodic(t[2][1],rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.],unitAngle='Degree')
test.testT(t,7)
# sur une base en radian (par defaut)
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 45., 5., (11,11,11))
C._initVars(a,'F',1.); C._initVars(a,'centers:G',2.)
C._addBC2Zone(a,'wall','BCWall','jmin')
C._addBC2Zone(a,'overlap','BCOverlap','jmax')
b = G.cylinder((0.,0.,0.), 0.1, 1., 45., 90., 5., (11,11,11)); b[0] = 'cyl2'
C._initVars(b,'F',1.); C._initVars(b,'centers:G',2.)
C._addBC2Zone(b,'wall','BCWall','jmin')
C._addBC2Zone(b,'overlap','BCOverlap','jmax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
C._addState(t[2][1], 'EquationDimension', 3)
t[2][1] = X.connectMatchPeriodic(t[2][1],rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,90.],unitAngle=None)
test.testT(t,8)

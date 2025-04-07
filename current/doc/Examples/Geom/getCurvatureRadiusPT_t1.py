# - getCurvatureRadius (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

test.TOLERANCE = 1.e-6

# cercle
a = D.circle((0,0,0), 1, 10, 0, 10)
a = D.getCurvatureRadius(a)
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
test.testT(t, 1)

# ligne : courbure infinie
a = D.line((0,0,0), (1,0,0), 3)
rad = D.getCurvatureRadius(a)
t = C.newPyTree(['Base',1]); t[2][1][2].append(rad)
test.testT(t,2)

# bezier
pts = D.polyline([(6,0.01,1), (5.4,0.036,1), (4.8,0.064,1), (2.5,0.21,1),
                  (0.3,0.26,1),(0,0.047,1),(0,0,0)])
a = D.bezier(pts, 100)
rad = D.getCurvatureRadius(a)
t = C.newPyTree(['Base',1]); t[2][1][2].append(rad)
test.testT([rad],3)

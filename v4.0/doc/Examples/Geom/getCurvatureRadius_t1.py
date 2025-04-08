# - getCurvatureRadius (array) -
import Geom as D
import Transform as T
import KCore.test as test

test.TOLERANCE = 1.e-6

# ligne : courbure infinie
a = D.line((0,0,0), (1,0,0), 3)
rad = D.getCurvatureRadius(a)
test.testA([rad],1)

# cercle
a = D.circle((0,0,0), 1, 10, 0, 10)
rad = D.getCurvatureRadius(a)
test.testA([rad],2)

# bezier
pts = D.polyline([(6,0.01,1), (5.4,0.036,1), (4.8,0.064,1), (2.5,0.21,1),
                  (0.3,0.26,1),(0,0.047,1),(0,0,0)])
a = D.bezier( pts, 100 )
rad = D.getCurvatureRadius(a)
test.testA([rad],3)

# cercle pas dans le plan
a = D.circle((0,0,0), 1, 10, 0, 10)
a = T.rotate(a, (0,0,0), (1,0,0), 32.)
rad = D.getCurvatureRadius(a)
test.testA([rad],4)

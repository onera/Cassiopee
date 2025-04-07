# - connect1D (array) -
import Geom as D
import KCore.test as test

# Test avec deux droites
P1 = [-0.5,0,0]; P1b = [0.5,0,0]
P2 = [1,-1.5,0]; P2b = [1,-0.5,0]
l1 = D.line(P1, P1b, N=10)
l2 = D.line(P2, P2b, N=10)

out = D.connect1D([l1,l2], sharpness=0)
test.testA([out], 1)

out = D.connect1D([l1,l2], sharpness=1)

# Test avec deux splines
p1 = D.polyline([(0,0,0), (1,1,0), (2,0,0)])
s1 = D.spline(p1, N=10)

p2 = D.polyline([(2.1,-0.1,0), (4,-2,0), (3,-3,0)])
s2 = D.spline(p2, N=10)
out = D.connect1D([s1,s2], sharpness=0, lengthFactor=0.5)
test.testA([out], 3)

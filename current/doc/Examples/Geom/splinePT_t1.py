# - spline (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# Spline 1D
c = D.polyline([(0.,0.,0.), (1.,1.,0.), (2.,1.,0.), \
                (3.,0.,0.), (4.,-1.,0.), (5.,6.,0.), \
                (6.,1.,0.), (7.,2.,0.), (8.,1.,0.), \
                (9.,-1.,0.), (10.,1.,0.), (11.,-1.,0.)])
d = D.spline(c,3,100); d[0] = 'spline'
t = C.newPyTree(['Base',1]); t[2][1][2].append(d)
test.testT(t, 1)

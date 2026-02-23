# - spline (array) -
import Generator as G
import Converter as C
import Geom as D

# Spline 1D
c = D.polyline([(0.,0.,0.), (1.,1.,0.), (2.,1.,0.), \
                (3.,0.,0.), (4.,-1.,0.), (5.,6.,0.), \
                (6.,1.,0.), (7.,2.,0.), (8.,1.,0.), \
                (9.,-1.,0.), (10.,1.,0.), (11.,-1.,0.)])
# With a specified number of points
d = D.spline(c, 3, N=100)
# With a specified density of points
e = D.spline(c, 3, density=10.)
C.convertArrays2File([c, d, e], 'out.plt')

# Spline 2D
ni = 4; nj = 4
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))

C.setValue(a, (1,1,1), [1.,1.,2.])
C.setValue(a, (1,2,1), [1.,2.,5.])
C.setValue(a, (1,3,1), [1.,3.,5.])
C.setValue(a, (1,4,1), [1.,4.,2.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,5.])
C.setValue(a, (2,4,1), [2.,4.,2.])
C.setValue(a, (3,1,1), [3.,1.,2.])
C.setValue(a, (3,2,1), [3.,2.,5.])
C.setValue(a, (3,3,1), [3.,3.,5.])
C.setValue(a, (3,4,1), [3.,4.,2.])
C.setValue(a, (4,1,1), [4.,1.,2.])
C.setValue(a, (4,2,1), [4.,2.,5.])
C.setValue(a, (4,3,1), [4.,3.,5.])
C.setValue(a, (4,4,1), [4.,4.,2.])

b = D.spline(a, 4, N=30, M=30)
c = D.spline(a, 4, density=10.)
C.convertArrays2File([a, b, c], 'out2.plt')

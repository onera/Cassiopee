# - bezier (array) -
import Geom as D
import Converter as C
import Generator as G

# Bezier 1D
pts = D.polyline([(0.,0.,0.), (0.,1.,0.), (2.,1.,0.), (2.,0.,0.),\
                  (4.,-1.,0.), (5.,6.,0.)])
# With a specified number of points
a = D.bezier(pts, N=100)
# With a specified point density
b = D.bezier(pts, density=10.)
C.convertArrays2File([pts, a, b], 'out.plt')

# Bezier 2D
ni = 2; nj = 3
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a, (1,1,1), [1.,1.,2.])
C.setValue(a, (1,2,1), [1.,2.,4.])
C.setValue(a, (1,3,1), [1.,3.,2.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,2.])
b = D.bezier(a, density=10.)
C.convertArrays2File([a]+[b], 'out2.plt')

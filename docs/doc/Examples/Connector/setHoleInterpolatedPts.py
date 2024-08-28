# - setHoleInterpolatedPts (array) -
import Converter as C
import Connector as X
import Generator as G

def sphere(x,y,z):
    if x*x+y*y+z*z < 0.5**2 : return 0.
    else: return 1.

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
celln = C.node2Center(a)
celln = C.initVars(celln, 'cellN', sphere, ['x','y','z'])
celln = X.setHoleInterpolatedPoints(celln, depth=1)
C.convertArrays2File([celln], 'out.plt')

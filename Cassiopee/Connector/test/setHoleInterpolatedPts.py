# - setHoleInterpolatedPts (array) -
import Converter as C
import Connector as X
import Generator as G

def sphere(x,y,z):
    return (x*x + y*y + z*z >= 0.5**2).astype(float)

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
celln = C.node2Center(a)
celln = C.initVars(celln, 'cellN', sphere, ['x','y','z'], isVectorized=True)
celln = X.setHoleInterpolatedPoints(celln, depth=1)
C.convertArrays2File([celln], 'out.plt')

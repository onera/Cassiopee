# - setInterpData (pyTree) -
# donneur NGON + leastsquare
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test

# NGON donor zone
a = G.cartNGon((0,0,-0.2),(0.01,0.01,0.1),(101,101,5))

# receiver
pts = D.circle((0.5,0.5,0),0.05,N=20)
C._initVars(pts, 'centers:cellN', 2)

pts2 = X.setInterpData(pts, a, order=3,
                       nature=1, loc='centers', storage='direct', hook=None,
                       method='leastsquares')
test.testT(pts2, 1)

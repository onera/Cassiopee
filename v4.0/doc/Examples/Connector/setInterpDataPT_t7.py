# - setInterpData (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test

# 2D structured donor zone
a = G.cart((0,0,0.),(0.01,0.01,0.1),(101,101,1))
pts = D.circle((0.5,0.5,0),0.05,N=20)
C._initVars(pts, 'cellN', 2); C._initVars(pts, '{centers:cellN}=2')
pts2 = X.setInterpData(pts, a, order=2,storage='direct',loc='nodes')
test.testT(pts2,1)
pts2 = X.setInterpData(pts, a, order=2,storage='direct',loc='centers')
test.testT(pts2,2)

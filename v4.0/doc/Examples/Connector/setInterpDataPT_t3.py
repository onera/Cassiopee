# - setInterpData (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test

# Tetra donor zone
a = G.cartTetra((0,0,-0.2),(0.01,0.01,0.1),(101,101,5))
points = D.circle((0.5,0.5,0),0.05,N=20)
points = C.convertArray2Tetra(points)

# Receiver 'BAR'
pts = C.initVars(points, 'cellN', 2) # interpole tout
pts2 = X.setInterpData(pts, a, loc='nodes')
test.testT(pts2,1)

pts = C.initVars(points, 'centers:cellN', 2) # interpole tout
pts2 = X.setInterpData(pts, a, loc='centers')
test.testT(pts2,2)

# Receiver 'NODE'
pts = C.initVars(points, 'cellN', 2) # interpole tout
pts = C.convertArray2Node(pts)
pts2 = X.setInterpData(pts, a, loc='nodes')
test.testT(pts2,3)

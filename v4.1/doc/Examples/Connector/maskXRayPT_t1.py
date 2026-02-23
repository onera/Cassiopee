# - maskXRay (pyTree) -
# sortie points de percage
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

# cas sphere 3D
surf = D.sphere((0,0,0), 0.5, 20)
surf = T.rotate(surf,(0.,0.,0.),(0.,1.,0.),90.)
surf = C.convertArray2Tetra(surf)
xray = X.maskXRay__(surf)
test.testT(xray,1)

# cas surface 2D avec body en BAR
surf = D.circle((0,0,0), 0.5, 0., 360.)
surf = C.convertArray2Tetra(surf)
xray =  X.maskXRay__(surf, 0.,2)
test.testT(xray,2)

# cas surface 2D avec body en TRI
surf = G.cylinder((0.,0.,0.), 0., 1., 360., 0., 1., (50,50,2))
surf = T.subzone(surf,(1,50,1),(50,50,2))
surf = C.convertArray2Tetra(surf)
xray =  X.maskXRay__(surf, 0.,2)
test.testT(xray,3)

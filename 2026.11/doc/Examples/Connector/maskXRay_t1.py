# - maskXRay (array) -
import KCore.test as test
import Connector as X
import Generator as G
import Geom as D
import Transform as T
import Converter as C

# cas sphere 3D
surf = D.sphere((0,0,0), 0.5, 20)
surf = T.rotate(surf,(0.,0.,0.),(0.,1.,0.),90.)
surf = C.convertArray2Tetra(surf)
res0 =  X.maskXRay__([surf])
test.testA([res0],1)

# cas surface 2D avec body en BAR
surf = D.circle((0,0,0), 0.5, 0., 360.)
surf = C.convertArray2Tetra(surf)
res0 =  X.maskXRay__([surf], 0.,2)
test.testA([res0],2)

# cas surface 2D avec body en TRI
surf = G.cylinder((0.,0.,0.), 0., 1., 360., 0., 1., (50,50,2))
surf = T.subzone(surf,(1,50,1),(50,50,2))
surf = C.convertArray2Tetra(surf)
res0 =  X.maskXRay__([surf], 0.,2)
test.testA([res0],3)

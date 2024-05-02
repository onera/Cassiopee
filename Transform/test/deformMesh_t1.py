# - deformMesh (array) -
import Transform as T
import Converter as C
import Geom as D
import KCore.test as test

a1 = D.sphere6((0,0,0), 1, 18)
a1 = C.convertArray2Tetra(a1); a1 = T.join(a1)
point = C.getValue(a1, 0)
a2 = T.deformPoint(a1, point, (0.1,0.05,0.2), 0.5, 2.)
delta = C.addVars(a1, ['dx','dy','dz'])
delta = C.extractVars(delta, ['dx','dy','dz'])
delta[1][:,:] = a2[1][:,:]-a1[1][:,:]
C._addVars([a1, delta])

test.stdTestA(T.deformMesh, a1)

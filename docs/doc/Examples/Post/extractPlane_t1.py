# - extractPlane (array) -
import Converter as C
import Generator as G
import Post as P
import Transform as T
import KCore.test as test

# test structure
m = G.cylinder((0.,0.,0.), 0., 1., 0, 360, 1., (50,20,3))
m = C.initVars(m, 'F', 1.)
m = T.subzone(m, (1,m[3],1), (m[2],m[3],m[4]))
m = T.reorder(m, (1,3,2))
m = G.close(m)
p = P.extractPlane([m], (0.,0.,1.,-0.5))
test.testA([p],1)

# test non structure tetra
m2 = C.convertArray2Tetra(m)
m2 = C.initVars(m2, 'F', 3.)
p = P.extractPlane([m2],(0.,0.,1.,-0.5))
test.testA([p],2)

# test non structure hexa
m2 = C.convertArray2Hexa(m)
p = P.extractPlane([m2],(0.,0.,1.,-0.5))
test.testA([p],3)

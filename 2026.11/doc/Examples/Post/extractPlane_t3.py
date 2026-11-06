# - extractPlane (array) -
import Converter as C
import Post as P
import Geom as D
import KCore.test as test

# test structure
m = D.sphere((0.,0.,0.),1.,N=20)
m = C.initVars(m, 'F', 1.)
p = P.extractPlane([m], (0.,0.,1.,-0.5))
test.testA([p],1)

# test non structure TRI
m2 = C.convertArray2Tetra(m)
m2 = C.initVars(m2, 'F', 3.)
p = P.extractPlane([m2],(0.,0.,1.,-0.5))
test.testA([p],2)

# test non structure QUAD
m2 = C.convertArray2Hexa(m)
p = P.extractPlane([m2],(0.,0.,1.,-0.5))
test.testA([p],3)

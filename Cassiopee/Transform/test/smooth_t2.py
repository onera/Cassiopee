# - smooth (array) -
import Transform as T
import Generator as G
import Post as P
import KCore.test as test

# Avec contraintes

# TRI
a = G.cartTetra((0.,0.,0.), (1,1,1), (10,10,1))
b = P.exteriorFaces(a)
a = T.smooth(a, eps=0.5, niter=10, fixedConstraints=b, delta=1.)
test.testA([a],1)

# QUAD
a = G.cartHexa((0.,0.,0.), (1,1,1), (10,10,1))
b = P.exteriorFaces(a)
a = T.smooth(a, eps=0.5, niter=10, fixedConstraints=b, delta=1.)
test.testA([a],2)

# STRUCT 2D
a = G.cart((0.,0.,0.), (1,1,1), (10,10,1) )
b = P.exteriorFaces(a)
a = T.smooth(a, eps=0.5, niter=10, fixedConstraints=b, delta=1.)
test.testA([a], 3)

# TETRA
a = G.cartTetra((0.,0.,0.), (1,1,1), (10,10,10) )
b = P.exteriorFaces(a)
a = T.smooth(a, eps=0.5, niter=10, fixedConstraints=b, delta=1.)
test.testA([a],4)

# HEXA
a = G.cartHexa((0.,0.,0.), (1,1,1), (10,10,10) )
b = P.exteriorFaces(a)
a = T.smooth(a, eps=0.5, niter=10, fixedConstraints=b, delta=1.)
test.testA([a],5)

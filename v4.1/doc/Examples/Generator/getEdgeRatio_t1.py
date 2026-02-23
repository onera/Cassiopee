# - getEdgeRatio (array) -
import Generator as G
import KCore.test as test

# STRUCT 2D
a = G.cart((0,0,0),(0.1,1,1),(11,11,1))
r = G.getEdgeRatio(a,dim=2)
test.testA([r],1)
# 2 plans en k
a = G.cart((0,0,0),(0.1,1,1),(11,11,2))
r = G.getEdgeRatio(a,dim=2)
test.testA([r],2)

# STRUCT 3D
a = G.cart((0,0,0),(0.1,1,1),(11,11,11))
r = G.getEdgeRatio(a)
test.testA([r],3)

# TRI
a = G.cartTetra((0,0,0),(0.1,1,1),(11,11,1))
r = G.getEdgeRatio(a)
test.testA([r],4)

# QUAD
a = G.cartHexa((0,0,0),(0.1,1,1),(11,11,1))
r = G.getEdgeRatio(a)
test.testA([r],5)

# TETRA
a = G.cartTetra((0,0,0),(0.1,1,1),(11,11,11))
r = G.getEdgeRatio(a)
test.testA([r],6)

# HEXA
a = G.cartHexa((0,0,0),(0.1,1,1),(11,11,11))
r = G.getEdgeRatio(a)
test.testA([r],7)

# PENTA
a = G.cartPenta((0,0,0),(0.1,1,1),(11,11,11))
r = G.getEdgeRatio(a)
test.testA([r],8)

# PYRA
a = G.cartPyra((0,0,0),(0.1,1,1),(11,11,11))
r = G.getEdgeRatio(a)
test.testA([r],9)

# NGON
a = G.cartNGon((0,0,0),(0.1,1,1),(11,11,11))
r = G.getEdgeRatio(a)
test.testA([r],10)

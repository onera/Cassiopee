# - addNormalLayers (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

# Tests sans lissage

d = C.array('d', 3, 1, 1);
d[1][0,0] = 0.1; d[1][0,1] = 0.2; d[1][0,2] = 0.3

# Structured (i,j-array)
a = D.sphere( (0,0,0), 1, 50 )
a = G.addNormalLayers(a, d)
test.testA([a], 1)

# Structured (i,j-array) avec champ
a = D.sphere( (0,0,0), 1, 50 )
a = C.initVars(a, 'Density', 1.)
a = G.addNormalLayers(a, d)
test.testA([a], 11)

# Unstructured (TRI)
a = D.sphere( (0,0,0), 1, 50 )
a = C.convertArray2Tetra(a)
a = G.addNormalLayers(a, d)
a = C.convertArray2Tetra(a)
test.testA([a], 2)

# Unstructured (TRI) avec champ
a = D.sphere( (0,0,0), 1, 50 )
a = C.convertArray2Tetra(a)
a = C.initVars(a, 'Density', 1.)
a = G.addNormalLayers(a, d)
a = C.convertArray2Tetra(a)
test.testA([a], 21)

# Unstructured (QUAD)
a = D.sphere( (0,0,0), 1, 50 )
a = C.convertArray2Hexa(a)
a = G.addNormalLayers(a, d)
a = C.convertArray2Tetra(a)
test.testA([a], 3)

# Structure (i-array)
a = D.line((0,0,0),(10,1,0))
a = G.addNormalLayers(a, d)
test.testA([a], 4)

#
# list of arrays
#
# Structured (i,j-array)
a = D.sphere( (0,0,0), 1, 50 )
res = G.addNormalLayers([a], d)
test.testA(res, 5)

# Unstructured (TRI)
a = D.sphere( (0,0,0), 1, 50 )
a = C.convertArray2Tetra(a)
res = G.addNormalLayers([a], d)
test.testA(res, 6)

# Unstructured (QUAD)
a = D.sphere( (0,0,0), 1, 50 )
a = C.convertArray2Hexa(a)
res = G.addNormalLayers([a], d)
test.testA(res, 7)

# Structure (i-array)
a = D.line((0,0,0),(10,1,0))
res = G.addNormalLayers([a], d)
test.testA(res, 8)

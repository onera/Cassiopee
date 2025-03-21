# - concatenate (PyTree) -

import Converter.PyTree   as C
import Generator.PyTree   as G
import Connector.PyTree   as X
import Transform.PyTree   as T
import Intersector.PyTree as XOR
import KCore.test         as test

# ----------------------------------------------------------------
# TEST 1
# ----------------------------------------------------------------
# Maillages
a    = G.cartHexa((0.,0.,0.), (1.,1.,1.), (3,2,2))
a    = C.convertArray2NGon(a)
a    = C.initVars(a,'{centers:varA} = 5.')
b    = G.cartHexa((2.,0.,0.), (1.,1.,1.), (3,2,2))
b    = C.convertArray2NGon(b)
b    = C.initVars(b,'{centers:varB} = 9.')

# BCs
subz = G.cartNGon( (0,0,0), (1,1,1), (3,2,1) )
hook = C.createHook(a, function='faceCenters')
ids  = C.identifyElements(hook, subz)
a    = C.addBC2Zone(a, 'wall', 'BCWall', faceList=ids)
subz = G.cartNGon( (2,0,0), (1,1,1), (3,2,1) )
hook = C.createHook(b, function='faceCenters')
ids  = C.identifyElements(hook, subz)
b    = C.addBC2Zone(b, 'wall', 'BCWall', faceList=ids)

# Raccord
t    = C.newPyTree(['Base',a,b])
t    = X.connectMatch(t)

# BCs
t    = C.fillEmptyBCWith(t, 'nref','BCFarfield')

# Champs
t   = C.initVars(t,'{centers:var} = 3.')

# Concatenation
t = XOR.concatenate(t)
test.testT(t, 1)

# ----------------------------------------------------------------
# TEST 2
# ----------------------------------------------------------------
# Maillages
a    = G.cartHexa((0.,0.,0.), (1.,1.,1.), (11,7,3))
a    = C.convertArray2NGon(a)

# BCs
subz = G.cartNGon( (0,0,0), (1,1,1), (11,7,1) )
hook = C.createHook(a, function='faceCenters')
ids  = C.identifyElements(hook, subz)
a    = C.addBC2Zone(a, 'wall', 'BCWall', faceList=ids)

t    = T.splitNParts(a, 6, multigrid=0, dirs=[1,2,3])

t    = C.initVars(t,'{centers:var} = 5.')

# Raccord
t    = X.connectMatch(t)

# BCs
t    = C.fillEmptyBCWith(t, 'nref','BCFarfield')

# Concatenation
t = XOR.concatenate(t)
test.testT(t, 2)

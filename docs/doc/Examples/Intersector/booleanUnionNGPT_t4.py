# - boolean union (PyTree) -
import Intersector.PyTree as XOR
import Connector.PyTree   as X
import Generator.PyTree   as G
import Converter.PyTree   as C
import Transform.PyTree   as T
import KCore.test         as test


h1 = G.cartHexa((0.,0.,0.), (1.,1.,1.), (13,12,7))
h2 = G.cartHexa((2.,1.,-0.5), (0.5,0.5,0.5), (23,17,33))

h1 = C.convertArray2NGon(h1); h1 = G.close(h1)
h2 = C.convertArray2NGon(h2); h2 = G.close(h2)

# BCs
subz = G.cartNGon( (0,0,0), (1,1,1), (13,12,1) )
hook = C.createHook(h1, function='faceCenters')
ids  = C.identifyElements(hook, subz)
h1   = C.addBC2Zone(h1, 'wall', 'BCWall', faceList=ids)

subz = G.cartNGon( (0,0,0), (1,1,1), (1,12,7) )
ids  = C.identifyElements(hook, subz)
h1   = C.addBC2Zone(h1, 'inflow', 'BCInflow', faceList=ids)

h1   = C.fillEmptyBCWith(h1, 'outflow','BCOutflow')

subz = G.cartNGon( (2.,1.,-0.5), (0.5,0.5,0.5), (23,17,1) )
hook = C.createHook(h2, function='faceCenters')
ids  = C.identifyElements(hook, subz)
h2   = C.addBC2Zone(h2, 'farfield', 'BCFarfield', faceList=ids)

subz = G.cartNGon( (2.,1.,-0.5), (0.5,0.5,0.5), (1,17,33) )
hook = C.createHook(h2, function='faceCenters')
ids  = C.identifyElements(hook, subz)
h2   = C.addBC2Zone(h2, 'sym', 'BCSymmetryPlane', faceList=ids)

h2   = C.fillEmptyBCWith(h2, 'far2','BCFarfield')

h1 = C.fillEmptyBCWith(h1, 'wall','BCWall')
h2 = C.fillEmptyBCWith(h2, 'nref','BCFarfield')

# Split
h1 = T.splitNParts(h1, 5, multigrid=0, dirs=[1,2,3])
h2 = T.splitNParts(h2, 4, multigrid=0, dirs=[1,2,3])

h1 = X.connectMatch(h1)
h2 = X.connectMatch(h2)

h1 = C.initVars(h1,'{centers:var1} = 3.')
h2 = C.initVars(h2,'{centers:var2} = 4.')
h1 = C.initVars(h1,'{centers:var}  = 5.')
h2 = C.initVars(h2,'{centers:var}  = 6.')

x = XOR.booleanUnion(h1, h2, tol=0, multi_zone=True)
test.testT(x,1)

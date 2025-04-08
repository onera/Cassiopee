# - writeZones (pyTree) -
# with Filter
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Filter as Filter
import KCore.test as test

LOCAL = test.getLocal()

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'centers:Density=1.')
b = G.cart((11,0,0), (1,1,1), (10,10,10))
C._initVars(b, 'centers:Density=1.')
t = C.newPyTree(['Base'])
C.convertPyTree2File(t, LOCAL+'/out.hdf')

h = Filter.Handle(LOCAL+'/out.hdf')

# Interface sur arbre (mais a doit etre mis dans t)
t[2][1][2] = [a,b]
h.writeZones(t, znp='Base/cart')

# Interface sur zones
h.writeZones(a, znp='Base/cart')
h.writeZones([a,b], znp=['Base/cart','Base/cart.0'])

r = C.convertFile2PyTree(LOCAL+'/out.hdf')
test.testT(r, 1)

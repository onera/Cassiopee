# - writeVariables (pyTree) -
# with Filter
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Filter as Filter

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, 'centers:Density=1.')
b = G.cart((11,0,0), (1,1,1), (10,10,10))
C._initVars(b, 'centers:Density=1.')
t = C.newPyTree(['Base'])
C.convertPyTree2File(t, 'out.hdf')

h = Filter.Handle('out.hdf')
h.writeZonesWoVars(a, znp='Base/cart')
h.writeVariables(a, 'centers:Density', znp='Base/cart')

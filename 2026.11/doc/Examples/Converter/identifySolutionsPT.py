# - identifySolutions (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

# A mesh with a field
ni = 21; nj = 21; nk = 21
m = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
m = C.initVars(m, '{Density}={CoordinateX}')

# Create extraction mesh
a = D.sphere((0,0,0),0.1)

# Identify solutions
hook = C.createGlobalHook([m], 'nodes')
a = C.identifySolutions(a, m, hookN=hook)
C.freeHook(hook)
C.convertPyTree2File(a, 'out.cgns')

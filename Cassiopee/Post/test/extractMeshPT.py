# - extractMesh (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(m, 'Density', 2.)
C._initVars(m, 'centers:cellN', 1)

# Extraction mesh
a = G.cart((0.,0.,0.5), (1., 0.1, 1.), (20, 20, 1)); a[0] = 'extraction'

# Extract solution on extraction mesh
a = P.extractMesh(m, a)
t = C.newPyTree(['Solution', 3, 'Extraction', 2])
t[2][1][2].append(m); t[2][2][2] += [a]
C.convertPyTree2File(a, 'out.cgns')

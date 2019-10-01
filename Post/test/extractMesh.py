# - extractMesh (array) -
import Converter as C
import Post as P
import Generator as G

ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', 1.)
# Create extraction mesh
a = G.cart((0.,0.,0.), (1., 0.1, 0.1), (20, 20, 1))
# Extract solution on extraction mesh
a2 = P.extractMesh([m], a)
C.convertArrays2File([m,a2], 'out.plt')

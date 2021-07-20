# - extractPoint (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

ni = 10; nj = 10; nk = 10;
a = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
def F(x,y,z): return x + 2.*y + 3.*z
a = C.initVars(a, 'F', F, ['CoordinateX','CoordinateY','CoordinateZ'])

# Utilisation directe
val = P.extractPoint(a, (0.55, 0.38, 0.12), 2); print(val)

# Utilisation avec un hook
hook = C.createHook(a, function='extractMesh')
val = P.extractPoint(a, (0.55, 0.38, 0.12), 2, hook=[hook]); print(val)

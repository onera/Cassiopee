# - extractPoint (array) -
import Converter as C
import Generator as G
import Post as P

ni = 10; nj = 10; nk = 10;
a = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
def F(x,y,z): return x*x*x*x + 2.*y + z*z
a = C.initVars(a, 'F', F, ['x','y','z'])

# Utilisation directe
val = P.extractPoint([a], (0.55, 0.38, 0.12), 2); print(val)

# Utilisation avec un hook
hook = C.createHook([a], function='extractMesh')
val = P.extractPoint([a], (0.55, 0.38, 0.12), 2, hook=hook); print(val)

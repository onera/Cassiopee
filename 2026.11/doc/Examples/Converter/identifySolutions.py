# - identifySolutions (array) -
import Converter as C
import Generator as G
import Geom as D
ni = 21; nj = 21; nk = 21
m = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
hook = C.createGlobalHook([m], function='nodes')
sol = C.initVars(m, 'ro={x}')
sol = C.extractVars(sol,['ro'])

# Create extraction mesh
a = D.sphere((0,0,0),0.1)
# Identify solutions of sol in a
a2 = C.identifySolutions(a, sol, hook)
C.freeHook(hook)
a = C.addVars([a,a2])
m = C.addVars([m,sol])
C.convertArrays2File([m,a], 'out.plt')

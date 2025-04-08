# - identifySolutions (array) -
import Converter as C
import Generator as G
import Geom as D
import KCore.test as test

# test: several variables
# Donor mesh structure
ni = 21; nj = 21; nk = 21
m = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
hook = C.createGlobalHook([m],function='nodes')
sol = C.initVars(m, '{rou}={x}')
sol = C.initVars(sol, '{rov}={y}')
sol = C.initVars(sol, '{row}={z}')

a = D.sphere((0,0,0),0.1)
a2 = C.identifySolutions(a,sol,hook,vars=[],tol=1000.)
test.testA([a2],1)
a2 = C.identifySolutions(a,sol,hook,vars=['rou'],tol=1000.)
test.testA([a2],2)

# - identifySolutions (array) -
import Converter as C
import Generator as G
import Geom as D
import Transform as T
import KCore.test as test

# Donor mesh structure
ni = 21; nj = 21; nk = 21
m = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
hook = C.createGlobalHook([m],function='nodes')
sol = C.initVars(m, '{ro}={x}')
sol = C.extractVars(sol,['ro'])

# RCV Structure
a = D.sphere((0,0,0),0.1)
a2 = C.identifySolutions(a,sol,hook,tol=1000.)
test.testA([a2],1)
# RCV TRI
a1 = C.convertArray2Tetra(a)
a2 = C.identifySolutions(a1,sol,hook,tol=1000.)
test.testA([a2],2)
# RCV NGON
a1 = C.convertArray2NGon(a)
a2 = C.identifySolutions(a1,sol,hook,tol=1000.)
test.testA([a2],3)
#
C.freeHook(hook)
#
# Dnrs: list of zones
#
ni = 21; nj = 21; nk = 21
m = T.splitSize(m,N=501)
nzones = len(m)
m[2] = C.convertArray2Tetra(m[2])
m[10] = C.convertArray2Hexa(m[10])
m[nzones-1] = C.convertArray2NGon(m[nzones-1])
hook = C.createGlobalHook(m,function='nodes')
sol = C.initVars(m, '{ro}={x}')
sol = C.extractVars(sol,['ro'])
a = D.sphere((0,0,0),0.1)
# RCV Structure
a2 = C.identifySolutions(a,sol,hook,tol=1000.)
test.testA([a2],4)
# RCV TRI
a1 = C.convertArray2Tetra(a)
a2 = C.identifySolutions(a1,sol,hook,tol=1000.)
test.testA([a2],5)
# NGON
a1 = C.convertArray2NGon(a)
a2 = C.identifySolutions(a1,sol,hook,tol=1000.)
test.testA([a2],6)
C.freeHook(hook)

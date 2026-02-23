# - enforceZ -
import Generator as G
import KCore.test as test

# Distribution : exact nb of added points
Ni = 50; Nj = 1; Nk = 50
a = G.cart((0,0,0), (1./(Ni-1), 1., 0.5/(Nk-1)), (Ni,Nj,Nk))
b = G.enforceZ(a, 0.3, 0.001, 10,15)
test.testA([b],1)

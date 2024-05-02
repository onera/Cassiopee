# - enforcePoint -
import Generator as G
import KCore.test as test

# distribution
Ni = 20; Nj = 20
a = G.cart((0,0,0), (1./(Ni-1),5./(Nj-1),1), (Ni,Nj,1))
b = G.enforcePoint(a, 0.5)
test.testA([b],1)

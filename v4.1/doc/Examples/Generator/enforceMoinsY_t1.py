# - enforceMoinsY (array) -
import Generator as G
import KCore.test as test

Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
# Monotonic distribution
b = G.enforceMoinsY(a, 1.e-3, (10,15))
test.testA([b],1)
# Exact nb of added pts
b = G.enforceMoinsY(a, 1.e-3, 10,15)
test.testA([b],2)

# - enforceMoinsX (array) -
import Generator as G
import KCore.test as test

# Distribution
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
# monotonic
b = G.enforceMoinsX(a, 1.e-3, (10,15))
test.testA([b],1)
# exact nb of added pts
b = G.enforceMoinsX(a, 1.e-3,10,15)
test.testA([b],2)

# - optimizeOverlap (array) -
import Converter as C
import Generator as G
import Transform as T
import Connector as X
import KCore.test as test

Ni = 50; Nj = 50; Nk = 2
a = G.cart((0,0,0),(1./(Ni-1), 1./(Nj-1),1), (Ni,Nj,Nk))
b = G.cart((0,0,0),(2./(Ni-1), 2./(Nj-1),1), (Ni,Nj,Nk))
a = T.rotate(a, (0,0,0), (0,0,1), 10.)
a = T.translate(a, (0.5,0.5,0))

ca = C.node2Center(a); ca =  C.initVars(ca, 'cellN', 1.)
cb = C.node2Center(b); cb =  C.initVars(cb, 'cellN', 1.)
# critere de volume
res = X.optimizeOverlap(a,ca,b,cb)
test.testA(res,1)
# priorites
res = X.optimizeOverlap(a, ca, b, cb, prio1=0, prio2=0)
test.testA(res,2)

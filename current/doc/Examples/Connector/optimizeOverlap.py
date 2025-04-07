# - optimizeOverlap (array) -
import Converter as C
import Generator as G
import Transform as T
import Connector as X

Ni = 50; Nj = 50; Nk = 2
a = G.cart((0,0,0),(1./(Ni-1), 1./(Nj-1),1), (Ni,Nj,Nk))
b = G.cart((0,0,0),(2./(Ni-1), 2./(Nj-1),1), (Ni,Nj,Nk))
a = T.rotate(a, (0,0,0), (0,0,1), 10.)
a = T.translate(a, (0.5,0.5,0))

ca = C.node2Center(a); ca =  C.initVars(ca, 'cellN', 1.)
cb = C.node2Center(b); cb =  C.initVars(cb, 'cellN', 1.)
res = X.optimizeOverlap(a, ca, b, cb)
C.convertArrays2File(res, "out.plt")

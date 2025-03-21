# - deform2 (array) -
import Converter as C
import Generator as G
import Transform as T
import time

# Cas 1
# a = G.cartNGon((0,0,0), (1,1,1), (5,5,1))
# C.convertArrays2File([a], 'out1.plt')
# #n = G.getNormalMap(a)
# #n = C.center2Node(n); n[1] = n[1]*10
# sx = C.initVars(a, '{sx}={x}*{x}')
# sy = C.initVars(a, '{sy}={y}')
# sz = C.initVars(a, 'sz', 1.)
# n = C.addVars([sx,sy,sz])
# n = C.extractVars(n, ['sx','sy','sz'])
# b = T.transform.deform2(a, n)
# print b
# b = G.close(b)
# C.convertArrays2File([a,b], 'out.plt')

# Morceau de sphere
#A = D.sphereYinYang((0,0,0), 1., N=20)
#a = A[0]

# Sphere a probleme
#a = D.sphere((0,0,0), 1., N=20)
#a = C.convertArray2Hexa(a)
#a = G.close(a)

# Sphere6 sans probleme
#a = D.sphere6((0,0,0), 1., N=20)
#a = C.convertArray2Hexa(a)
#a = T.join(a); a = G.close(a)

# Le cube
#a = G.cart((0,0,0), (1,1,1), (10,10,10))
#a = C.convertArray2Hexa(a)
#a = P.exteriorFaces(a)

# Le toit
a = G.cart((0,0,0), (1,1,1), (10,10,1))
b = T.rotate(a, (0,0,0), (0,1,0), 120.)
b = T.translate(b, (9,0,0))
c = [a,b]
c = C.convertArray2Hexa(c)
c = T.join(c)
a = G.close(c)
#C.convertArrays2File(c, 'out1.plt')

# go
a = C.convertArray2NGon(a)
a = G.close(a)
a = T.reorder(a, (+1,))
ainit = C.copy(a)

out = [ainit]
for i in range(4):
    n = G.getNormalMap(a)
    n = C.center2Node(n)
    n = C.normalize(n, ['sx','sy','sz'])
    n[1] = n[1]*0.3
    b = T.transform.deform2(a, n)
    b = G.close(b, 1.e-6)
    out.append(b)
    #b = T.reorder(b, (1,))

    #CPlot.display([ainit,b])
    a = C.copy(b)
    time.sleep(0.1)
    if a[2].size > 1e6: break

C.convertArrays2File(out, 'out.plt')

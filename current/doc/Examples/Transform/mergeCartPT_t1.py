# - mergeCart (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

dh = 0.1; n = 11
A = []
a1 = G.cart((0.,0.,0.),(dh,dh,dh),(n,n,n)); a1 = T.oneovern(a1,(2,2,2))
a2 = G.cart((1.,0.,0.),(dh,dh,dh),(n,n,n)); a2 = T.oneovern(a2,(2,2,2))
a3 = G.cart((1.,1.,0.),(dh,dh,dh),(n,n,n))
a4 = G.cart((0.,1.,0.),(dh,dh,dh),(n,n,n))
A = [a1,a2,a3,a4]
for i in range(1,5): A[i-1][0] = 'cart'+str(i)

t = C.newPyTree(['Base']); t[2][1][2] += A
t[2][1][2] = T.mergeCart(t[2][1][2])
test.testT(t)

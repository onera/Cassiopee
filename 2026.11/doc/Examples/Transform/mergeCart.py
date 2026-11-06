# - mergeCart (array) -
import Converter as C
import Generator as G
import Transform as T

dh = 0.1; n = 11
A = []
a1 = G.cart((0.,0.,0.),(dh,dh,dh),(n,n,n)); A.append(a1)
a2 = G.cart((1.,0.,0.),(dh,dh,dh),(n,n,n)); A.append(a2)
a3 = G.cart((1.,1.,0.),(dh,dh,dh),(n,n,n)); A.append(a3)
a4 = G.cart((0.,1.,0.),(dh,dh,dh),(n,n,n)); A.append(a4)
A[0] = T.oneovern(A[0],(2,2,2))
A[1] = T.oneovern(A[1],(2,2,2))
res = T.mergeCart(A)
C.convertArrays2File(res, "out.plt")

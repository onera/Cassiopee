# - integNormProduct -
import Converter as C
import Post as P
import Transform as T
from numpy import *

# Lit le fichier et le met dans arrays
arrays = []

ar = C.convertFile2Arrays("cyl01.tp", "fmt_tp")
arrays.append(ar[0])
ar = C.convertFile2Arrays("cyl02.tp", "fmt_tp")
arrays.append(ar[0])
ar = C.convertFile2Arrays("cyl03.tp", "fmt_tp")
arrays.append(ar[0])
ar = C.convertFile2Arrays("cyl04.tp", "fmt_tp")
arrays.append(ar[0])
ar = C.convertFile2Arrays("cyl05.tp", "fmt_tp")
arrays.append(ar[0])
ar = C.convertFile2Arrays("cyl06.tp", "fmt_tp")
arrays.append(ar[0])

n = len(arrays)
a = range(0, n)
data = []
A=[]
for a in arrays :
    ni = a[2]
    nj = a[3]
    nk = a[4]
    print(ni, nj, nk)
    dens = ones( (1, (ni-1) * (nj-1) * (nk-1)), float64 )
#    dens = ones( (3, (ni-1) * (nj-1) * (nk-1)), float64 )
    densa = ['t', dens, (ni-1), (nj-1), (nk-1)]  # en centres
    ac = T.subzone(a,(1,1,1),(ni,nj,1))
    data.append(densa)
    A.append(ac)

res = P.integ(A,data,[])
#res = P.integNormProduct(A, data,[])
print(res)

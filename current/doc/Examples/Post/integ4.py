# - integ
# - integNorm
import Converter as C
import Post as P
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

#C.convertFile(arrays,"mesh6Blocks.plt","bin_tp")
n = len(arrays)
a = range(0, n)
data = []
arraysc=[]

for i in a :
    ni = arrays[i][2]
    nj = arrays[i][3]
    nk = arrays[i][4]
    print(ni, nj, nk)
#    dens = ones( (1, ni * nj * nk), float64 )
#    densa = ['t', dens, ni, nj, nk]
    dens = ones( (1, (ni-1) * (nj-1) * (nk-1)), float64 )
    densa = ['t', dens, (ni-1), (nj-1), (nk-1)]
    data.append(densa)
    c=P.node2Center(arrays[i])
    arraysc.append(c)

print(arraysc[0])
print(data[0])
res = P.integ([arraysc[0]], [data[0]], [])
#print(res)
#res = P.integNorm(data)
print(res)

# - integNormProduct -
import Converter as C
import Transform as T
import Post as P
from numpy import *

# Lit le fichier et le met dans arrays
arrays = []

ar = C.convertFile2Arrays("cyl01.tp", "fmt_tp")
a2 = T.subzone(ar[0], (1,1,1), (101,301,1))
arrays.append(a2)
ar = C.convertFile2Arrays("cyl02.tp", "fmt_tp")
a2 = T.subzone(ar[0], (1,1,1), (51,301,1))
arrays.append(a2)
ar = C.convertFile2Arrays("cyl03.tp", "fmt_tp")
a2 = T.subzone(ar[0], (1,1,1), (71,91,1))
arrays.append(a2)
ar = C.convertFile2Arrays("cyl04.tp", "fmt_tp")
a2 = T.subzone(ar[0], (1,1,1), (71,91,1))
arrays.append(a2)
ar = C.convertFile2Arrays("cyl05.tp", "fmt_tp")
a2 = T.subzone(ar[0], (1,1,1), (491,211,1))
arrays.append(a2)
ar = C.convertFile2Arrays("cyl06.tp", "fmt_tp")
a2 = T.subzone(ar[0], (1,1,1), (491,211,1))
arrays.append(a2)

n = len(arrays)
a = range(0, n)
density = []

for i in a :
    ni = arrays[i][2]
    nj = arrays[i][3]
    nk = arrays[i][4]
    dens = ones( (3, (ni-1) * (nj-1) * nk), float64 )
    densa = ['t', dens, (ni-1), (nj-1), nk]
    density.append(densa)

res = P.integNormProduct(arrays, density, [])
print(res)
res = P.integMoment(arrays, density, [], (0.,0.,0.))
print(res)

# - integNorm (array) -
import Generator as G
import Converter as C
import Transform as T
import Post as P
from numpy import *
import math

res1 = math.pi*10.
res1c = (math.pi-math.pi/25.)*10.
res2 = 2.
res2c = 1.62
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,2,30))
a2 = T.subzone(a, (1,1,1), (50,1,30))

# integNorm node2center, nj = 1
ni = a[2]-1; nj = a[3]-1; nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]
res = P.integNorm([a2],[densa],[])
res = res[0]
if math.fabs(res[0]) > 1.e-1:
    print("pb in integNormNodeCenter, nj=1")

# integNorm, nj = 1
ni = a[2]; nj = a[3]-1; nk = a[4]
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]
res = P.integNorm([a2],[densa],[])
res = res[0]
if math.fabs(res[0]) > 1.e-1:
    print("pb in integNorm, nj=1")

##############################################################
a = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
a2 = T.subzone(a, (1,1,1), (11,11,1))
#C.convertArrays2File([a],'out.plt','bin_tp')

# integNorm node2center, nk = 1

ni = a[2]-1
nj = a[3]-1
nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]
res = P.integNorm([a2],[densa],[])
res = res[0]
if math.fabs(res[2]-res2) > 1.e-1:
    print("pb in integNormNodeCenter, nk=1")

# integNorm, nk = 1

ni = a[2]
nj = a[3]
nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]

res = P.integNorm([a2],[densa],[])
res = res[0]
if math.fabs(res[2]-res2) > 1.e-1:
    print("pb in integNorm, nk=1")


#######################################################
a = G.cart( (0,0,0), (1., 0.2, 0.1), (2, 11, 11))
a2 = T.subzone(a, (1,1,1), (1,11,11))
#C.convertArrays2File([a],'out.plt','bin_tp')

# integNorm node2center, ni = 1
ni = a[2]-1
nj = a[3]-1
nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]

res = P.integNorm([a2],[densa],[])
res = res[0]
if math.fabs(res[0]-res2) > 1.e-1:
    print("pb in integNormNodeCenter, ni=1")


# integ, ni = 1
ni = a[2]-1
nj = a[3]
nk = a[4]
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]
res = P.integNorm([a2],[densa],[])
res = res[0]
if math.fabs(res[0]-res2) > 1.e-1:
    print("pb in integ, ni=1")

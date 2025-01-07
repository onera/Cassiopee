# - integ -
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

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (100,2,30))
a2 = T.subzone(a, (1,1,1), (100,1,30))
#C.convertArrays2File([a], "new.plt", "bin_tp")

# integ node2center, nj = 1
ni = a[2]-1; nj = a[3]-1; nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]
res = P.integ([a2],[densa],[])
if math.fabs(res[0]-res1) > 1.e-1:
    print("pb in integNodeCenter, nj=1")

# integ, nj = 1
ni = a[2]; nj = a[3]-1; nk = a[4]
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]
res = P.integ([a2],[densa],[])
if math.fabs(res[0]-res1) > 1.e-1:
    print("pb in integ, nj=1")

##############################################################
a = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
a2 = T.subzone(a, (1,1,1), (11,11,1))
#C.convertArrays2File([a],'out.plt','bin_tp')

# integ node2center, nk = 1

ni = a[2]-1; nj = a[3]-1; nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]

res = P.integ([a2],[densa],[])
if math.fabs(res[0]-res2) > 1.e-1:
    print("pb in integNodeCenter, nk=1")

# integ, nk = 1

ni = a[2]; nj = a[3]; nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]

res = P.integ([a2],[densa],[])
if math.fabs(res[0]-res2) > 1.e-1:
    print("pb in integ, nk=1")

#######################################################
a = G.cart( (0,0,0), (1., 0.2, 0.1), (2, 11, 11))
a2 = T.subzone(a, (1,1,1), (1,11,11))
#C.convertArrays2File([a],'out.plt','bin_tp')

# integ node2center, ni = 1
ni = a[2]-1
nj = a[3]-1
nk = a[4]-1
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]

res = P.integ([a2],[densa],[])
if math.fabs(res[0]-res2) > 1.e-1:
    print("pb in integNodeCenter, ni=1")


# integ, ni = 1
ni = a[2]-1; nj = a[3]; nk = a[4]
dens = ones( (1, ni*nj*nk), float64 )
densa = ['t', dens, ni,nj,nk]

res = P.integ([a2],[densa],[])
if math.fabs(res[0]-res2) > 1.e-1:
    print("pb in integ, ni=1")

print('done.')

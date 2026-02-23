# - integNormProduct (array) -
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
# Lit le fichier et le met dans arrays
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,2,30))
a2 = T.subzone(a, (1,1,1), (50,1,30))
#C.convertArrays2File([a], "new.plt", "bin_tp")

# integNormProd node2center, nj = 1
ni = a[2]-1; nj = a[3]-1; nk = a[4]-1
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integNormProduct([a2],[densa],[])
if math.fabs(res) > 1.e-1:
    print("pb in integNormProdNodeCenter, nj=1")

# integNormProd, nj = 1
ni = a[2]
nj = a[3]-1
nk = a[4]
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integNormProduct([a2],[densa],[])
if math.fabs(res) > 1.e-1:
    print("pb in integNormProd, nj=1")

##############################################################
a = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
a2 = T.subzone(a, (1,1,1), (11,11,1))
#C.convertArrays2File([a],'out.plt','bin_tp')

# integNormProd node2center, nk = 1

ni = a[2]-1
nj = a[3]-1
nk = a[4]-1
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integNormProduct([a2],[densa],[])
if math.fabs(res) > 1.e-1:
    print("pb in integNormProdNodeCenter, nk=1")

# integNormProd, nk = 1
ni = a[2]
nj = a[3]
nk = a[4]-1
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integNormProduct([a2],[densa],[])
if math.fabs(res) > 1.e-1:
    print("pb in integNormProd, nk=1")

#######################################################
a = G.cart( (0,0,0), (1., 0.2, 0.1), (2, 11, 11))
a2 = T.subzone(a, (1,1,1), (1,11,11))
#C.convertArrays2File([a],'out.plt','bin_tp')

# integNormProd node2center, ni = 1
ni = a[2]-1
nj = a[3]-1
nk = a[4]-1
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integNormProduct([a2],[densa],[])
if math.fabs(res) > 1.e-1:
    print("pb in integNormProdNodeCenter, ni=1")


# integNormProd, ni = 1
ni = a[2]-1
nj = a[3]
nk = a[4]
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integNormProduct([a2],[densa],[])
if math.fabs(res) > 1.e-1:
    print("pb in integNormProd, ni=1")
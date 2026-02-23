# - integMoment (array) -
import Generator as G
import Converter as C
import Transform as T
import Post as P
from numpy import *
import math
#
res1 = math.pi*10.
res1c = (math.pi-math.pi/25.)*10.
res2 = 1.
res2c = 1.62
#
#------------------------------------------------------------
#
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,2,30))
a2 = T.subzone(a, (1,1,1), (50,1,30))
C.convertArrays2File([a,a2], "new.plt", "bin_tp")
#
# integMoment node2center, nj = 1
#
ni = a2[2]-1; nj = a2[3]; nk = a2[4]-1
dens = C.array('tx,ty,tz', ni, nj, nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integMoment([a2],[densa],[], (0.,0.,5.))
print(res, res1)


# integMoment, nj = 1
ni = a2[2]
nj = a2[3]
nk = a2[4]
dens = C.array('tx,ty,tz',ni,nj,nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integMoment([a2],[densa],[], (0.,0.,5.))
## if math.fabs(res[0]-res1) > 1.e-1:
##     print "pb in integMoment, nj=1"
print(res, res1)

##############################################################
a = G.cart( (0,0,0), (0.1, 0.1, 1.), (11, 11, 2))
a2 = T.subzone(a, (1,1,1), (11,11,1))
C.convertArrays2File([a2],'out.plt','bin_tp')

# integMoment node2center, nk = 1

ni = a2[2]
nj = a2[3]
nk = a2[4]
dens = C.array('tx,ty,tz',ni,nj,nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integMoment([a2],[densa],[], (0.5,1.,0.))
## if math.fabs(res[0]-res2) > 1.e-1:
##     print "pb in integMomentNodeCenter, nk=1"
print(res, res2)

# integMoment, nk = 1

ni = a[2]
nj = a[3]
nk = a[4]-1
dens = C.array('tx,ty,tz',ni,nj,nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)
res = P.integMoment([a2],[densa],[], (0.5,1.,0.))
print(res, res2)

#######################################################
a = G.cart( (0,0,0), (1., 0.2, 0.1), (2, 11, 11))
a2 = T.subzone(a, (1,1,1), (1,11,11))
C.convertArrays2File([a2],'out.plt','bin_tp')

# integMoment node2center, ni = 1
ni = a[2]-1
nj = a[3]-1
nk = a[4]-1
dens = C.array('tx,ty,tz',ni,nj,nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)

res = P.integMoment([a2],[densa],[], (0.,1.,0.5))
## if math.fabs(res[0]-res2) > 1.e-1:
##     print "pb in integMomentNodeCenter, ni=1"
print(res, res2)

# integMoment, ni = 1
ni = a2[2]; nj = a2[3]; nk = a2[4]
dens = C.array('tx,ty,tz',ni,nj,nk)
densa = C.initVars(dens,'tx', 0.)
densa = C.initVars(densa,'ty', 1.)
densa = C.initVars(densa,'tz', 0.)

res = P.integMoment([a2],[densa],[], (0.,1.,0.5))
## if math.fabs(res[0]-res2) > 1.e-1:
##     print "pb in integMoment, ni=1"
print(res, res2)

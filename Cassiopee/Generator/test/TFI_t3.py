# - TFI tri
# - TFI penta
import Generator as G
import Geom as D
import KCore.test as test

#------------
# 1- TFI tri
#------------
l1 = D.line((0,0,0),(0,1,0), 15)
l2 = D.line((0,0,0),(1,0,0), 15)
l3 = D.line((1,0,0),(0,1,0), 15)
tri = G.TFI([l1,l2,l3])
test.testA([tri],1)

#--------------
# 2- TFI penta
#--------------
n = 5
P1 = (0,0,0)
P2 = (1,0,0)
P3 = (0,1,0)

P4 = (0,0,1)
P5 = (1,0,1)
P6 = (0,1,1)

# Tmin
l1 = D.line(P1,P3, n)
l2 = D.line(P1,P2, n)
l3 = D.line(P2,P3, n)
tri1 = G.TFI([l1,l2,l3])

# Tmax
l1 = D.line(P4,P6, n)
l2 = D.line(P4,P5, n)
l3 = D.line(P5,P6, n)
tri2 = G.TFI([l1,l2,l3])

# imin
p = 30
l1 = D.line(P1,P3, n)
l2 = D.line(P4,P6, n)
l3 = D.line(P1,P4, p)
l4 = D.line(P3,P6, p)
quad1 = G.TFI([l3,l4,l1,l2])

# jmin
l1 = D.line(P1,P2, n)
l2 = D.line(P4,P5, n)
l3 = D.line(P1,P4, p)
l4 = D.line(P2,P5, p)
quad2 = G.TFI([l3,l4,l1,l2])

# diag
l1 = D.line(P2,P3, n)
l2 = D.line(P5,P6, n)
l3 = D.line(P2,P5, p)
l4 = D.line(P3,P6, p)
quad3 = G.TFI([l3,l4,l1,l2])
penta = G.TFI([tri1,tri2,quad1,quad2,quad3])
test.testA([penta],2)

#
# TFI TETRA
#
n = 15
P1 = (0,0,0); P2 = (1,0,0); P3 = (0,1,0); P4 = (0,0,1)
# face P1P2P3
l1 = D.line(P1,P2, n) #P1P2
l2 = D.line(P1,P3, n) #P1P3
l3 = D.line(P2,P3, n) #P2P3
tri1 = G.TFI([l2,l1,l3])

# face P2P3P4
l1 = D.line(P2,P3, n) #P2P3
l2 = D.line(P2,P4, n) #P2P4
l3 = D.line(P3,P4, n) #P3P4
tri2 = G.TFI([l2,l1,l3])

# face P3P4P1
l1 = D.line(P3,P4, n) #P3P4
l2 = D.line(P3,P1, n) #P3P1
l3 = D.line(P4,P1, n) #P4P1
tri3 = G.TFI([l2,l1,l3])

# face P4P1P2
l1 = D.line(P4,P1, n) #P4P1
l2 = D.line(P4,P2, n) #P4P2
l3 = D.line(P1,P2, n) #P1P2
tri4 = G.TFI([l2,l1,l3])

TRI = [tri1,tri2,tri3,tri4]
tetra = G.TFI(TRI)
test.testA([tetra],3)

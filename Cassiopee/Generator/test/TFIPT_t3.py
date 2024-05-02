# test des fonctions (pyTree)
# - TFI tri
# - TFI penta
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

#------------
# 1- TFI tri
#------------
l1 = D.line((0,0,0),(0,1,0), 15)
l2 = D.line((0,0,0),(1,0,0), 15)
l3 = D.line((1,0,0),(0,1,0), 15)
tri = G.TFI([l1,l2,l3])
test.testT(tri,1)

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

test.testT(penta,2)

# - TFI PENTA (array) -
import Generator as G
import Converter as C
import Geom as D

n = 5
P1 = (0,0,0); P2 = (1,0,0); P3 = (0,1,0)
P4 = (0,0,1); P5 = (1,0,1); P6 = (0,1,1)

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

out = [tri1,tri2,quad1,quad2,quad3]
C.convertArrays2File(out, 'out1.plt')

r = [tri1,tri2,quad1,quad2,quad3]
penta = G.TFI(r)

C.convertArrays2File(out+[penta], 'out2.plt')

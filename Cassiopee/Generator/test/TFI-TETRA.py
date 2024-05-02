# - TFI TETRA (array) -
import Generator as G
import Converter as C
import Geom as D
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
C.convertArrays2File([tetra],"out.plt")

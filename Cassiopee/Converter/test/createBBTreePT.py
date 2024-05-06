# - createBBTree (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter

b = G.cart((10,0,0), (1,1,1), (10,2,2))

t = C.newPyTree(['Base'])
for i in range(200):
    a = G.cart((i,0,0), (1,1,1), (2,2,2))
    t[2][1][2].append(a)

bbt = Cmpi.createBBTree(t)
ret = Cmpi.intersect(b, bbt); print(ret)
#>> [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
Converter.converter.deleteBBTree(bbt)

# - splitTBranches (array)
import Converter as C
import Generator as G
import Transform as T

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,50))
c1 = T.subzone(a,(1,1,1),(50,1,1))
c2 = T.subzone(a,(1,1,50),(50,1,50))
c3 = T.subzone(a,(1,1,1),(1,1,50))
c = [c1,c2,c3]; c = C.convertArray2Hexa(c)
c = T.join(c)
res = T.splitTBranches(c)
C.convertArrays2File(res,"out.plt")

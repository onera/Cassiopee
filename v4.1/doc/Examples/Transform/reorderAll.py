# - reorderAll (array) -
import Converter as C
import Generator as G
import Transform as T

ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = T.rotate(m1, (0.2,0.2,0.), (0.,0.,1.), 15.)
m2 = T.reorder(m2,(-1,2,3))
a = [m1,m2]
a = T.reorderAll(a,1)
C.convertArrays2File(a, "out.plt")

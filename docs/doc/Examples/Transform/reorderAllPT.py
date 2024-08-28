# - reorderAll (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

ni = 30; nj = 40; nk = 1
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk)); m1[0]='cart1'
m2 = T.rotate(m1, (0.2,0.2,0.), (0.,0.,1.), 15.)
m2 = T.reorder(m2,(-1,2,3)); m2[0]='cart2'
t = C.newPyTree(['Base',2,[m1,m2]])
t = T.reorderAll(t, 1)
C.convertPyTree2File(t, "out.cgns")

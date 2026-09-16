# - getCollidingCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Post.PyTree as P

t1 = G.cart((0,0,0), (1,1,1), (10,10,10))
t1 = C.convertArray2NGon(t1)
t2 = G.cart((1.,1.5,3.), (1,1,1), (10,10,10))
t2 = C.convertArray2NGon(t2)

res = XOR.getCollidingCells(t1, t2, RTOL=0.05)

[ids_in1,ids_in2] = res[0]
m = XOR.getCells(t1, [ids_in1], are_face_ids=False)

C.convertPyTree2File(m, "out.cgns")
C.convertPyTree2File(t1, "t1.cgns")
C.convertPyTree2File(t2, "t2.cgns")

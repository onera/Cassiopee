import Geom.PyTree as D
import Converter.PyTree as C
import Transform.PyTree as T
import Intersector.PyTree as XOR
import KCore.test as test

a = D.sphere6((0,0,0), 1, N=20)
b = a[0]

a = C.convertArray2NGon(a)
a = T.join(a)

t = C.newPyTree(['Base1',a,'Base2',b])

t1 = XOR.XcellN(t, [(0,1)], output_type=0)
test.testT(t1,1)
t2 = XOR.XcellN(t, [(0,1)], output_type=1)
test.testT(t2,2)
t3 = XOR.XcellN(t, [(0,1)], output_type=2)
test.testT(t3,3)
#C.convertPyTree2File(t3, "Sph611_0.cgns")

t1 = XOR.XcellN(t, [(1,0)], output_type=0)
test.testT(t1,4)
t2 = XOR.XcellN(t, [(1,0)], output_type=1)
test.testT(t2,5)
t3 = XOR.XcellN(t, [(1,0)], output_type=2)
test.testT(t3,6)
#C.convertPyTree2File(t3, "Sph612_0.cgns")

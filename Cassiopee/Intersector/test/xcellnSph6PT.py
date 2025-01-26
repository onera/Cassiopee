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
t2 = XOR.XcellN(t, [(0,1)], output_type=1)
t3 = XOR.XcellN(t, [(0,1)], output_type=2)

C.convertPyTree2File(t1, 'out1_0.cgns')
C.convertPyTree2File(t2, 'out2_0.cgns')
C.convertPyTree2File(t3, 'out3_0.cgns')

t1 = XOR.XcellN(t, [(1,0)], output_type=0)
t2 = XOR.XcellN(t, [(1,0)], output_type=1)
t3 = XOR.XcellN(t, [(1,0)], output_type=2)

C.convertPyTree2File(t1, 'out4_0.cgns')
C.convertPyTree2File(t2, 'out5_0.cgns')
C.convertPyTree2File(t3, 'out6_0.cgns')

# - sphereYinYang (PyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Converter.Internal as I
import Transform.PyTree as T
import time

s = D.sphereYinYang((0,0,0), 1., 50)
zs = I.getZones(s)
t = C.newPyTree(['Base', zs[0], 'Base2', zs[1]])

priorities = []
priorities.append((0,1))

t0 = time.time()
XOR._XcellN(t, priorities, output_type=2)
t1 = time.time()
print(' - XCELLN CPU time : ',t1-t0,'s')

C.convertPyTree2File(t, "out.cgns")


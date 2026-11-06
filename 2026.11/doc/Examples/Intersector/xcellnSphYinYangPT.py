# - XcellN (PyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Converter.Internal as I

s = D.sphereYinYang((0,0,0), 1., 50)
zs = I.getZones(s)
t = C.newPyTree(['Base', zs[0], 'Base2', zs[1]])

priorities = []
priorities.append((0,1))

XOR._XcellN(t, priorities, output_type=2)

C.convertPyTree2File(t, "out.cgns")

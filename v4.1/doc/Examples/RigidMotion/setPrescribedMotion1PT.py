# - setPrescribedMotion1 (pyTree) -
# Motion defined by time string
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion1(a, 'trans', tx="{t}")

C.convertPyTree2File(a, 'out.cgns')

# - evalPosition (pyTree) -
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion1(a, 'trans', tx="{t}")
b = R.evalPosition(a, time=0.1)

C.convertPyTree2File(b, 'out.cgns')

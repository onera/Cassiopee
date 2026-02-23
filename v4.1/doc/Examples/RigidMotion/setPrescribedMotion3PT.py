# - setPrescribedMotion3 (pyTree) -
# Motion defined by a constant rotation and translation speed
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion3(a, 'mot', transl_speed=(1,0,0))

C.convertPyTree2File(a, 'out.cgns')

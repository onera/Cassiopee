# - setPrescribedMotion2 (pyTree) - 
# Motion defined by a Cassiopee Solver rotor motion
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion2(a, 'rotor', transl_speed=(0.1,0,0), rot_omg=1.)

C.convertPyTree2File(a, 'out.cgns')

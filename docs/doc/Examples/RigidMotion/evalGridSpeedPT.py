# - evalGridSpeed (pyTree) -
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion3(a, 'motion', transl_speed=(1,0,0))
b = R.evalPosition(a, time=0.1)
R._evalGridSpeed(b, time=0.1)
C.convertPyTree2File(b, 'out.cgns')

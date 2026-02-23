# - evalPosition pour motion 3 (pyTree) -
# Motion defined by a constant rotation and translation speed
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion3(a, 'mot', transl_speed=(0.1,0,0),
                           axis_pnt=(1.,2.,0.), axis_vct=(0.,0.,1.),
                           omega=0.2)

b1 = R.evalPosition(a, time=1.); b1[0]='moved1'
b2 = R.evalPosition(a, time=2.); b2[0]='moved2'
b3 = R.evalPosition(a, time=3.); b3[0]='moved3'
C.convertPyTree2File([b1,b2,b3], 'out.cgns')

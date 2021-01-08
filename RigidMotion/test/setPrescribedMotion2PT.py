# - setPrescribedMotion2 (pyTree) - 
# Motion defined by a rotor motion
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Generator.PyTree as G

# Mime une pale suivant x, quart avant
a = G.cart((0.2,-0.075,0), (0.01,0.01,0.1), (131,11,1))
# Mettre tous les parametres
R._setPrescribedMotion2(a, 'rotor', transl_speed=(0.1,0,0), rot_omg=1.)
C.convertPyTree2File(a, 'out.cgns')

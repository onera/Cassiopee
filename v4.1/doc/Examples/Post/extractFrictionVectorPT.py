# - extractFrictionVector (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.ExtraVariables2 as PE

a = G.cart((0.,0.,0.), (13./100.,13./100.,1.), (100,100,1))
for n in ['ViscosityMolecular', 'Pressure',
          'gradxVelocityX', 'gradxVelocityY','gradxVelocityZ',
          'gradyVelocityX','gradyVelocityY','gradyVelocityZ',
          'gradzVelocityX','gradzVelocityY','gradzVelocityZ']:
    C._initVars(a, '{centers:%s} = 1.'%n)
PE._extractShearStress(a)
PE._extractFrictionVector(a)
C.convertPyTree2File(a, 'out.cgns')

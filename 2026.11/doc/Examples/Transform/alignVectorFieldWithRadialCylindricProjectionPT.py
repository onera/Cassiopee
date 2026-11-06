# - alignVectorFieldWithRadialCylindricProjection (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

surface = G.cart((-2.,-2.,0),(1,1,1),(5,5,1))
C._initVars(surface,'vx',0.)
C._initVars(surface,'vy',0.)
C._initVars(surface,'vz',1.)
T._alignVectorFieldWithRadialCylindricProjection(surface, (0,0,-1), (0,1,0), ['vx','vy','vz'])
C.convertPyTree2File(surface, 'out.cgns')

import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere6( (0,0,0), 1., N=20, ntype='QUAD')
C._initVars(a, 'Fx={CoordinateX}*0.1')
C._initVars(a, 'Fy={CoordinateY}*0.1')
C._initVars(a, 'Fz={CoordinateZ}*0.1')
C.convertPyTree2File(a, 'tryVector.cgns')

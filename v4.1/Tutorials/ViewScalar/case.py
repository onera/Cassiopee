import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere6( (0,0,0), 1., N=20, ntype='QUAD')
C._initVars(a, 'F={CoordinateX}')
C.convertPyTree2File(a, 'tryScalar.cgns')

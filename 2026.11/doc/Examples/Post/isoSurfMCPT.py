# - isoSurfMC (pyTree) -
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((-20,-20,-20), (0.5,0.5,0.5), (50,50,50))
a = C.initVars(a, '{field}={CoordinateX}*{CoordinateX}+{CoordinateY}*{CoordinateY}+{CoordinateZ}')

iso = P.isoSurfMC(a, 'field', value=5.)
C.convertPyTree2File(iso, 'out.cgns')

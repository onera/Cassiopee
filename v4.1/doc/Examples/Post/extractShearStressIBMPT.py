# - extractShearStress for IBM (pyTree) -
import Post.IBM as P_IBM
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Generator.PyTree as G

a = D.sphere((0,0,0),0.5,N=30)
a = C.convertArray2Tetra(a); a = G.close(a)
C._initVars(a,'{utau}={CoordinateX}**2')
C._initVars(a,'{VelocityX}={CoordinateZ}*{CoordinateY}')
C._initVars(a,'{VelocityY}={CoordinateX}*{CoordinateZ}')
C._initVars(a,'{VelocityZ}={CoordinateX}*{CoordinateY}')
C._initVars(a,'{Density}=1')

P_IBM._extractShearStress(a)
C.convertPyTree2File(a, "out.cgns")

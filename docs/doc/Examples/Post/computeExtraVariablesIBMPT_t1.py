#compute shear stress for IBM
import Post.IBM as P_IBM
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

a=D.sphere((0,0,0),0.5,N=30)
a = C.convertArray2Tetra(a); a = G.close(a)
C._initVars(a,'{centers:utau}={centers:CoordinateX}**2')
C._initVars(a,'{centers:VelocityX}={centers:CoordinateZ}*{centers:CoordinateY}')
C._initVars(a,'{centers:VelocityY}={centers:CoordinateX}*{centers:CoordinateZ}')
C._initVars(a,'{centers:VelocityZ}={centers:CoordinateX}*{centers:CoordinateY}')
C._initVars(a,'{centers:Density}=1')
C._initVars(a,'{centers:Pressure}=0.71')

P_IBM._computeExtraVariables(a,PInf=0.71, QInf=0.005, variables=['Cp','Cf','ShearStress'])
test.testT(a)

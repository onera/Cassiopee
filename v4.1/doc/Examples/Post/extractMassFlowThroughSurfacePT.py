import Converter.PyTree as C
import Post.IBM as P_IBM
import Geom.IBM as D_IBM
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal

a = D.cylinder((0.,0.,0.),1., 5, N=20,ntype='TRI')
a = T.splitSharpEdges(a)
tb = C.newPyTree(['Body'])
C._addFamily2Base(tb[2][1],'inlet')

for z in a:
    if C.getMaxValue(z,'CoordinateZ')==0:
        Internal.newFamilyName(name='FamilyName', value='inlet', parent=z)
tb[2][1][2]=a

h = 0.1
a = G.cart((-1.5,-1.5,-1.5),(h,h,h),(int(3/h)+1,int(3/h)+1,int(6.5/h)+1))
C._initVars(a,'centers:Density',1.)
C._initVars(a,'centers:VelocityX',0.)
C._initVars(a,'centers:VelocityY',0.)
C._initVars(a,'centers:VelocityZ',0.2)

massflow, tmf = P_IBM.extractMassFlowThroughSurface(tb, a, famZones=['inlet'])
print(massflow)
C.convertPyTree2File(tmf,"out.cgns")

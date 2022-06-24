# - setIBCType (pyTree) -
import Converter.PyTree as C
import Geom.IBM as D_IBM
import Geom.PyTree as D
import Transform.PyTree as T
import Converter.Internal as Internal

#Prepare Base
t=Internal.newCGNSBase('Base')
C._addFamily2Base(t, 'UP')
C._addFamily2Base(t, 'DOWN')

#Prepare Geometry
a = D.circle((0,0,0), 1. , 0., 360.)
a = D_IBM.setIBCType(a, "Musker")
a = T.splitNParts(a, 2, dirs=[1])
C._tagWithFamily(Internal.getNodeByName(a,'circle.1'), 'UP')
C._tagWithFamily(Internal.getNodeByName(a,'circle.3'), 'DOWN')
Internal.addChild(t, a)

#Set Snear based on familyName
t = D_IBM.setIBCType(t,"Log",familyName=['UP'])
C.convertPyTree2File(t, 'out.cgns')

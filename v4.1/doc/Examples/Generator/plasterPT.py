# - plaster (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Post.PyTree as P
import Geom.PyTree as D

a = D.sphere((0,0,0), 1)
a = T.subzone(a, (6,1,1), (50,200,1))
a = C.convertArray2Hexa(a); a = G.close(a)

# contours
c = P.exteriorFaces(a)
cs = T.splitConnexity(c)

# plaster hole
p = G.plaster([cs[0]],[a])
t = C.newPyTree(['Sphere',2,'Contours',1,'Plaster',2])
t[2][1][2].append(a); t[2][2][2] += cs; t[2][3][2].append(p)
C.convertPyTree2File(t, 'out.cgns')

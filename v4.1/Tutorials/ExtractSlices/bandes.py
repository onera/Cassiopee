# Extraction de bandes x=cte sur un maillage surfacique
import Converter.PyTree as C
import Post.PyTree as P
import Geom.PyTree as D
import Generator.PyTree as G
import Intersector.PyTree as Intersector
import Transform.PyTree as T

# IN: surface+solution
# Attention: les normales doivent etre dirigees vers l'exterieur
a = D.sphere( (0,0,0), 1, N=100)
a = C.initVars(a, 'p={CoordinateX}+{CoordinateY}')
a = C.convertArray2Tetra(a) ; a = G.close(a)

# Extraction de lignes x = cste
bb = G.bbox(a)
xmin = bb[0]; ymin = bb[1]; zmin = bb[2]
xmax = bb[3]; ymax = bb[4]; zmax = bb[5]

bands = []
dx = 1./10
for i in range(10):
        x = i/10.
        b = G.cart( (x,ymin-1,zmin-1), (dx,ymax-ymin+2,zmax-zmin+2), (2,2,2) )
        b = P.exteriorFaces(b) ; b = C.convertArray2Tetra(b)
        b = T.join(b) ; b = G.close(b)
        b = T.reorder(b, (-1,))
        band = Intersector.booleanIntersection(a, b)
        bands += [band]

# Boolean operators do not keep field
bands = P.extractMesh(a, bands)

# Sortie
t = C.newPyTree(['Base']) ; t[2][1][2] += bands
C.convertPyTree2File(t, 'out.cgns')

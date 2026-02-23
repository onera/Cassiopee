# Extraction de lignes x=cte sur un maillage surfacique
import Converter.PyTree as C
import Post.PyTree as P
import Geom.PyTree as D

# IN: surface+solution
a = D.sphere((0,0,0), 1, N=100)
a = C.initVars(a, 'p={CoordinateX}+{CoordinateY}')

# Extraction de lignes x = cste
lines = []
for i in range(10):
        x = i/10.
        p = P.isoSurfMC(a, 'CoordinateX', x)
        p[0][0] = "line%s"%str(i)
        lines += p

# Sortie
t = C.newPyTree(['sphere',a,
                 'lines', lines])
C.convertPyTree2File(t, 'out.cgns')

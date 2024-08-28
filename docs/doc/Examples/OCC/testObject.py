# test CAD objects
import OCC.PyTree as OCC
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

# cree et lit la CAD
cad = OCC.CAD('cube.step', 'fmt_step')

# Liste des objets faces et edges de la CAD
faces = cad.faces
edges = cad.edges

# test eval an edge
N = 10; h = 1./(N-1.)
d = G.cart((0,0,0), (h,1,1), (N,1,1))
m = edges[0].valueAt(d)
#C.convertPyTree2File(m, 'out.cgns')

# test eval an edge at one parameter
m = edges[0].valueAt(0.5)
#C.convertPyTree2File(m, 'out.cgns')

# test eval a face
d = G.cart((0,0,0), (h,h,1), (N,N,1))
m = faces[0].valueAt(d)
#C.convertPyTree2File(m, 'out.cgns')

# test eval face in one parameter point
m = faces[0].valueAt((0.5,0.5))
#C.convertPyTree2File(m, 'out.cgns')

# test de projection sur toutes les faces
z = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
cad._project(z)

# test de projection sur une liste de faces
cad._project(z, [1,2,3,4,5,6])

# test de projection de point
pt = D.point((20,-5,26))
cad._project(pt)
#C.convertPyTree2File(pt, 'out.cgns')

# test de projection de ligne
line = D.line((20,-5,26),(-10, 20, 0))
cad._project(line)
#C.convertPyTree2File(line, 'out.cgns')

# test des mailleurs, mesh all faces
m = cad.mesh('STRUCT', N=11)
C.convertPyTree2File(m, 'struct.cgns')

# test de l'association maillage/CAD
for z in m:
    face = cad.getLinkFace(z); print(face.number)

# test des autres mailleurs
m = cad.mesh('TRI', N=11)
C.convertPyTree2File(m, 'tri.cgns')

m = cad.mesh('QUADHO', N=11)
C.convertPyTree2File(m, 'quadho.cgns')

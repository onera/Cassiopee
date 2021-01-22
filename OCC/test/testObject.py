import OCC.PyTree as OCC
import Converter.PyTree as C

# cree et lit la CAD
cad = OCC.CAD('cube.step', 'fmt_step')

# Liste des objets faces et edges de la CAD
faces = cad.faces
edges = cad.edges

# test evalEdge
import Generator.PyTree as G
N = 10; h = 1./(N-1.)
d = G.cart((0,0,0), (h,1,1), (N,1,1))
m = cad.evalEdge(edges[0], d)
C.convertPyTree2File(m, 'out.cgns')

# test evalFace
d = G.cart((0,0,0), (h,h,1), (N,N,1))
m = cad.evalFace(faces[0], d)
C.convertPyTree2File(m, 'out.cgns')

# test de projection sur toutes les faces
z = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
cad._project(z)

# test de projection sur une liste de faces
cad._project(z, [1,2,3,4,5,6])

# test des mailleurs, mesh all faces
m = cad.mesh('STRUCT', N=11)
C.convertPyTree2File(m, 'struct.cgns')

# test de l'association maillage/CAD
for z in m:
    face = cad.getLinkFace(z); print(face.no)

# test des autres mailleurs
m = cad.mesh('TRI', N=11)
C.convertPyTree2File(m, 'tri.cgns')
m = cad.mesh('QUADHO', N=11)
C.convertPyTree2File(m, 'quadho.cgns')

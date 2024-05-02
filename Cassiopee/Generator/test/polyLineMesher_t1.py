# - polyLineMesher (array) -
import Generator.PolyLine as GP
import Generator as G
import Geom as D
import KCore.test as test
import Converter as C

a = [D.polyline([(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(0.,0.,0.)])]

m  = []
for i in a:
    i = C.convertArray2Tetra(i)
    i = G.close(i, 1.e-2)
    m.append(i)

# Donnees
h = 0.05; hf = 0.0001; density = 200

# Par familles
coords = []; walls = []; overlaps = []
for i in m:
    b = GP.polyLineMesher(i, h, hf, density)
    coords.append(b[0]); walls.append(b[1]); overlaps.append(b[2])

# A plat
meshes = []
for i in coords: meshes = meshes + i
test.testA(meshes,1)

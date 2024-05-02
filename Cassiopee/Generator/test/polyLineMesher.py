# - polyLineMesher (array) -
import Converter as C
import Generator.PolyLine as GP
import Generator as G
import Transform as T

# Read a 2D geometry created with tecplot
a = C.convertFile2Arrays('fusee.plt')
a = G.close(a,1e-2); a = T.reorder(a,(-1,2,3))
C.convertArrays2File(a, 'input.plt')

# Data
h = 0.02; hf = 0.0001; density = 500

# Per families
coords = []; walls = []
for i in a:
    b = GP.polyLineMesher(i, h, hf, density)
    coords.append(b[0])
    walls.append(b[1])

# Flat
meshes = []
for i in coords: meshes = meshes + i

C.convertArrays2File(meshes, 'out.plt')

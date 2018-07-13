# - indiceFace2Connect (array) -
import KCore
import Geom as D
import Converter as C
import Transform as T
import Generator as G
import numpy

# Maillage non structure
a = G.cartTetra((0,0,0), (1.,1,1), (10,10,1) )
a = G.close(a)

# liste des faces
elt = 0
faces = [elt*3+0, elt*3+1, elt*3+2]
connect = KCore.indiceFace2Connect(a, faces)
print connect
C.convertArrays2File([a], 'out.plt')

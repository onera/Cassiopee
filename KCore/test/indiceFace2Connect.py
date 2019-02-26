# - indiceFace2Connect (array) -
import KCore
import Converter as C
import Generator as G

# Maillage non structure
a = G.cartTetra((0,0,0), (1.,1,1), (10,10,1) )
a = G.close(a)

# liste des faces
elt = 0
faces = [elt*3+0, elt*3+1, elt*3+2]
connect = KCore.indiceFace2Connect(a, faces)
C.convertArrays2File([a], 'out.plt')

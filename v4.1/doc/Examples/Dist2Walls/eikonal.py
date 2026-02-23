# - eikonal (array) -
import Dist2Walls
import Generator as G
import Converter as C
import Geom as D
import Connector as X

# Bloc cartesien
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(128,128,1))
cellN = C.initVars(a, 'cellN', 1)
cellN = C.extractVars(cellN, ['cellN'])

# Init wall
sphere = D.sphere((6.4,6.4,0), 1., 30)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
cellN = X.blankCellsTri([a], [cellN], [sphere], blankingType=0)

a = C.addVars([a, cellN[0]])
a = C.initVars(a, '{Phi}=1.e12*({cellN}>0.)')
a = C.initVars(a,'{speed}=%f'%(1./0.1))
# Eikonal
b = Dist2Walls.eikonal(a)
C.convertArrays2File([b], 'out.plt')

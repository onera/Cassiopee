# Refine an octree mesh inside a surface defined by a sphere
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Connector.PyTree as X

# IN: sphere
a = D.sphere((0,0,0), 2, N=20)

# Octree mesh generation around the sphere
o = G.octree([a], [0.5], dfar=15., balancing=1)
t = C.newPyTree(['Base', o])
import numpy as N; BM = N.array([[1]])
def F(c, vol):
    if c == 0. and vol > volmin: return 1.
    else: return 0.

# Refine the octree inside the sphere until all the elements inside the
# sphere are of finest level
end = 0
while end == 0:
    t = X.blankCells(t, [[a]], BM, blankingType='center_in')
    t = G.getVolumeMap(t); volmin = C.getMinValue(t, 'centers:vol')
    t = C.initVars(t, 'centers:indicator', F, ['centers:cellN','centers:vol'])
    end = 1
    if  C.getMaxValue(t,'centers:indicator') == 1.: end = 0
    # Maintien du niveau de raffinement le plus fin
    o = t[2][1][2][0]; o = G.adaptOctree(o, 'centers:indicator', balancing=1)
    # Sortie
    t = C.newPyTree(['Base', o])

C.convertPyTree2File(t, 'out.cgns')

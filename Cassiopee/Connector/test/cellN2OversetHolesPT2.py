# - cellN2OversetHoles (pyTree) -
# - Dumping the OversetHoles node to files -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter

a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((0.5,0.5,0.5),(1,1,1),(10,10,10))
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)
t = C.addBC2Zone(t, 'overlap1', 'BCOverlap', 'imin')
C._initVars(t, 'centers:cellN', 0)
t = X.applyBCOverlaps(t)
t = X.cellN2OversetHoles(t)

zones = Internal.getNodesFromType(t, 'Zone_t')
for z in zones:
    ho = Internal.getNodesFromType(z, 'OversetHoles_t')
    if ho != []:
        h = ho[0][2][1][1]
        array = ['cell_index', h, h.size, 1, 1]
        Converter.convertArrays2File([array], 'hole_'+z[0]+'.v3d',
                                     'bin_v3d')

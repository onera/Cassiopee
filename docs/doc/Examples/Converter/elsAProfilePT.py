# - elsAProfile (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Converter.elsAProfile as CE
import Connector.PyTree as X
import Transform.PyTree as T
import numpy

#-----------------------------------------------------------------------------
# PARAMETRES
R0 = 1. # rayon du cylindre
R1 = 3.
ni = 201
DEPTH = 0 # Grilles recouvrantes ou non ?
vmin = 11 # nb de pts minimum par grille cartesienne
dhov = 1.e-1 # taille de maille pres des raccords chimere grilles de corps
#-----------------------------------------------------------------------------
a = G.cylinder((0,0,0),R0,R1,360,0,1,(ni,10,1)); a[0] = 'cylindre1'
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'ov1', 'BCOverlap', 'jmax')
#
t = C.newPyTree(['Cylindre1',a])
t = X.connectMatch(t,dim=2)

#---------------------------
# Octree mesh generation
#---------------------------
# frontieres exterieures des grilles de corps : surfaces de reference pour l ocrtree
extFaces = C.extractBCOfType(t,'BCOverlap')
# Taille de maille de l octree pres des surfaces de reference
dho = dhov*(vmin-1)
snears = [dho]*len(extFaces)
# generation de l octree
o = G.octree(extFaces, snears, dfar=10., balancing=1)
o = G.expandLayer(o)
# conversion en grilles cartesiennes patch grids
if DEPTH == 0: ext = 0
else: ext = DEPTH+1 # si Chimere
res = G.octree2Struct(o, vmin=11, ext=ext, optimized=1, merged=1)
res = C.fillEmptyBCWith(res,'nref','BCFarfield',dim=2)
#
# Ajout de la base cartesienne
t = C.addBase2PyTree(t,'CARTESIAN'); t[2][2][2] = res
# ajout d un plan en k
t = T.addkplane(t)
#
#---------------------------
# Assemblage Chimere
#---------------------------
DEPTH = 2
# bodies description
bodies = []
bases = Internal.getNodesFromType(t,'CGNSBase_t')
for b in bases:
    walls = C.extractBCOfType(b, 'BCWall')
    if walls != []: bodies.append(walls)

# blanking
BM = numpy.array([[0],[1]])
t = X.blankCells(t, bodies, BM,depth=DEPTH, dim=2)
t = X.setHoleInterpolatedPoints(t,depth=DEPTH)
t = X.applyBCOverlaps(t,depth=DEPTH)
t = X.optimizeOverlap(t,priorities=['Cylindre1',0,'CARTESIAN',1])
t = X.maximizeBlankedCells(t,depth=DEPTH)
t = X.setInterpolations(t,loc='cell')



# Arbre a la Cassiopee
C.convertPyTree2File(t, 'out_cassiopee.cgns')

# Arbre a la elsA
t = CE.convert2elsAxdt(t)
C.convertPyTree2File(t, 'out_elsa.cgns')
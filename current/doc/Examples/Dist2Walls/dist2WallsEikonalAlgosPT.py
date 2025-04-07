# - dist2WallsEikonal (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import Generator.PyTree as G
import numpy
import time

alg = DTW.fim

print("Preparation de test")
beg = time.time()
DEPTH = 2
# Bloc cartesien
N = 128; h = 0.1
# N = 4096; h = 0.025/8.
a = G.cart((0.,0.,0.),(h,h,h),(N,N,1))
# Init wall
sphere = D.sphere((N*h/2.,N*h/2.,0), 1., 100)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base',a])
C._initVars(t,'cellN', 1)
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='node_in')
# Condition aux limites
t = X.setHoleInterpolatedPoints(t,depth=1,loc='nodes')
C._initVars(t,'{flag}=({cellN}>1.)')
end = time.time()
print("Temps preparation test : {} secondes".format(end-beg))
print("Fin preparation du test. Appel solveur Eikonal")
beg =time.time()
t = DTW.distance2WallsEikonal(t,sphere,DEPTH=DEPTH,nitmax=10,algo=alg)
end = time.time()
print("Temps passe par python pour solveur Eikonal : {} secondes".format(end-beg))
C.convertPyTree2File(t, 'out.cgns')

# Bloc cartesien N=64
N = 64; h = 0.2
a = G.cart((0.,0.,0.),(h,h,h),(N,N,1)); a[0] = 'cart2'
# Init wall
#sphere = D.sphere((6.4,6.4,0), 1., 100)
sphere = D.sphere((N*h/2.,N*h/2.,0.), 1., 100)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
C._initVars(t,'cellN', 1.)
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='node_in')
t = X.setHoleInterpolatedPoints(t,depth=1,loc='nodes')
# Initialise le front
C._initVars(t,'{flag}=({cellN}>1.)')
t = DTW.distance2WallsEikonal(t,sphere,DEPTH=DEPTH,nitmax=10,algo=alg)
C.convertPyTree2File(t, 'out2.cgns')

# - gencartmb (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

# body mesh
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
h = 1.e-1  # Step of finest Cartesian grid
Dfar = 20. # Distance to far boundaries

# Nb of points per level:
# Here 4 levels, but last one is computed automatically
nlvl = [10,10,5] # nlvl[0]: coarse grid

t = C.newPyTree(['Bodies', 'CARTESIAN']); t[2][1][2].append(a)
zones = G.gencartmb(t[2][1], h, Dfar, nlvl)
t[2][2][2] += zones
C.convertPyTree2File(t, 'out.cgns')

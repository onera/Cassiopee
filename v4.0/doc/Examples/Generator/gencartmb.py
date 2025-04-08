# - gencartmb (array) -
import Generator as G
import Converter as C

# body grid
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
h = 1.e-1# Step of finest Cartesian grid
Dfar = 10.# Extension of far boundaries
nlvl = [5,5,5] # Nb of points per level, except the 4th level (automatic)
cartGrids = G.gencartmb([a], h, Dfar, nlvl)
C.convertArrays2File(cartGrids, 'out.plt')

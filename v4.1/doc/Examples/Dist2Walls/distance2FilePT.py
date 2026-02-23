# - distance2Walls (pyTree) -
# - Dump TurbulentDistance node to a file -
import Dist2Walls.PyTree as DW
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Internal as Internal
import Converter
import numpy

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
a = C.initVars(a,'centers:cellnf',1.)
sphere = D.sphere((1.2,0.,0.),0.2,100)
sphere = C.initVars(sphere,'centers:cellnf',1.)

t = C.newPyTree(['Base',a])
bodies = C.newPyTree(['Bodies',sphere])
t = DW.distance2Walls(t, bodies)
nodes = Internal.getNodesFromName(t, 'TurbulentDistance')
c = 0
for n in nodes:
    ni = n[1].shape[0]; nj = n[1].shape[1]; nk = n[1].shape[2]
    a = numpy.reshape(n[1], (ni*nj*nk), order='Fortran')
    a = numpy.reshape(a, (1,ni*nj*nk))
    array = ['walldistance', a, ni, nj, nk]
    array = Converter.initVars(array, 'wallglobalindex', 1)
    Converter.convertArrays2File([array], 'dist'+str(c)+'.v3d',
                                 'bin_v3d')
    c += 1

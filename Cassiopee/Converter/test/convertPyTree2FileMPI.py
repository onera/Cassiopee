# - convertPyTree2File (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Generator.PyTree as G

if Cmpi.rank == 0:
    a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
    a[0] = 'cart0'
else:
    a = G.cart((10,0,0),(0.1,0.1,0.1),(11,11,11))
    a[0] = 'cart1'

t = C.newPyTree(['Base', a])
Cmpi.convertPyTree2File(t, 'out.cgns')

# - generateIBMMesh (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.IBM as G_IBM
import Geom.IBM as D_IBM

tb = D.sphere((0.,0.,0.), 1., N=100)
tb = D_IBM.setDfar(tb, 10.)
tb = D_IBM.setSnear(tb, 0.02)
tb = C.newPyTree(['Base', tb])

t = G_IBM.generateIBMMesh(tb, dimPb=3, vmin=21, octreeMode=1, check=False)

C.convertPyTree2File(t, 'out.cgns')

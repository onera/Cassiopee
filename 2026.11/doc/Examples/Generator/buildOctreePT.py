# - buildOctree (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.IBM as G_IBM

tb = D.sphere((0.,0.,0.), 1., N=100)
tb = C.newPyTree(['Base', tb])

to = G_IBM.buildOctree(tb, dimPb=3, vmin=21, snears=0.02, dfars=10., tbox=None, octreeMode=1)

C.convertPyTree2File(to, 'out.cgns')
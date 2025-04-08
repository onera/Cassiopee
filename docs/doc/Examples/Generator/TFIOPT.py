# - TFIO (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

a = D.circle((0,0,0), 1., N=41)
r = G.TFIO(a)
C.convertPyTree2File(r, 'out.cgns')

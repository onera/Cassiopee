# - setIBCType (pyTree) -
import Converter.PyTree as C
import Geom.IBM as D_IBM
import Geom.PyTree as D

a = D.circle((0,0,0), 1. , 0., 360.)
a = D_IBM.setIBCType(a,"Musker")

C.convertPyTree2File(a, 'out.cgns')
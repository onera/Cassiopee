# - nurbs (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G

ni = 10; nj = 10
a = G.cart((0,0,0), (1,1,1), (ni,nj,1));
C._initVars(a,'weight',1.)
C.setValue(a,'weight',(7,1,1), 7.)
C.setValue(a,'weight',(9,5,1), 9.)
d = D.nurbs(a,'weight',4,100,100)
C.convertPyTree2File(d, 'out.cgns')

a = D.polyline ([(4.1,0.1,1.1),(1.1,0.2,1.2),(1.1,1.3,1.3),(1.1,1.5,1.4),(4.5,2.5,1.5),(5.6,1.5,1.6),(6.7,1.7,1.7),(7.8,0.8,1.8),(8.9,-1.9,1.9),(9,0,1)])
a = C.initVars(a,'weight',1.)
C.setValue(a, 'weight', (7,1,1), 7.)
C.setValue(a, 'weight', (9,1,1), 9.)
b = D.nurbs(a,'weight',4,2000)
C.convertPyTree2File(b, 'out2.cgns')

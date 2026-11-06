# - projectDir (array) -
import Geom as D
import Converter as C
import Generator as G
import Transform as T
a = D.sphere((0,0,0), 1., 20)
b = G.cart((1.1,-0.1,-0.1),(0.03,0.03,0.03), (1,50,50))
c = T.projectDir(b, [a], (1.,0,0))
d = T.projectDir([b], [a], (1.,0,0), smooth=1)
C.convertArrays2File([a,b,c]+d, 'out.plt')

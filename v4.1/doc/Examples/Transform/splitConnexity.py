# - splitConnexity (array) -
import Converter as C
import Transform as T
import Geom as D

a = D.text2D("CASSIOPEE")
B = T.splitConnexity(a)
C.convertArrays2File(B, 'out.plt')

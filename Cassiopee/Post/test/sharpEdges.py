# - sharpEdges ( array) -
import Converter as C
import Generator as G
import Post as P
import Transform as T
a1 = G.cart((0.,0.,0.),(1.5,1.,1.),(2,2,1))
a2 = T.rotate(a1,(0.,0.,0.),(0.,1.,0.),100.)
res = P.sharpEdges([a1,a2],alphaRef=45.)
C.convertArrays2File([a1,a2]+res,"out.plt")

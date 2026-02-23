# - convertArrays2File (fmt_svg) -
import Converter as C
import Generator as G

a = C.convertFile2Arrays('test.svg', 'fmt_svg', nptsCurve=130)
C.convertArrays2File(a, 'out.plt', 'bin_tp')
C.convertArrays2File(a, 'out.svg', 'fmt_svg')

a = G.cart( (0,0,0), (1.,1.,1.), (10, 10, 1) )
C.convertArrays2File([a], 'out0.svg', 'fmt_svg')

a = G.cart( (0,0,0), (100,100,100), (3, 3, 3) )
b = G.cart( (400,0,-2), (100,100,1), (3, 3, 1) )
c = G.cart( (0,300,0),(100,100,100), (3,3,3) )
d = G.cart( (400,300,-2),(100,100,1), (3,3,1) )
import Transform as T
a = T.rotate(a, (0,0,0), (1,1,0), 20.)
c = T.rotate(c, (0,300,0), (1,1,0), 20.)
C.convertArrays2File([a,b,c,d], 'out1.svg', 'fmt_svg')

b = C.convertArray2Tetra(b)
d = C.convertArray2Hexa(d)
C.convertArrays2File([a,b,c,d], 'out2.svg', 'fmt_svg')

a = C.convertArray2Tetra(a)
c = C.convertArray2Hexa(c)
C.convertArrays2File([a,b,c,d], 'out3.svg', 'fmt_svg')

# - cylinder (array) -
import Generator as G
import Converter as C
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
C.convertArrays2File([a], "out.plt")

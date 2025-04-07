# Generate input files for Converter tests
# Require Geom, Generator, Converter modules
import Converter as C
import Generator as G

cyl = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
cart = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
C.convertArrays2File([cyl, cart], "in.plt", "bin_tp")

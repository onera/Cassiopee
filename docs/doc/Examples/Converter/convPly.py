# - convertFile2Arrays (binary PLY) -
import Converter as C
import Generator as G

a = G.cartHexa( (0,0,0), (1,1,1), (10,10,1) )
C.convertArrays2File([a], 'out.ply')
a = C.convertFile2Arrays('out.ply')

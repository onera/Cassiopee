# - convertFile2Arrays (binary python pickle) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (100,100,100))
C.convertArrays2File([a], 'out.pickle', 'bin_pickle')
b = C.convertFile2Arrays('out.pickle', 'bin_pickle')

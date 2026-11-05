# - convertFile2Arrays (fmt free) -
import Converter as C
import Generator as G

a = C.array('x,y,f', 10, 1, 1)
a[1][0,:] = 1.
a[1][1,:] = 2.
C.convertArrays2File(a, 'out.csv')

a = C.convertFile2Arrays('out.csv')
print(a)

# - convertFile2Arrays (binaire tecplot) -
# - convertArrays2File (binaire tecplot) -
# Extrait certaines variables d'un fichier
# Ecrit un nouveau fichier tecplot
import Converter as C
import numpy as N

# Read a file into numpy arrays
arrays = C.convertFile2Arrays('in.plt')

out = []

for b in arrays:
    t = b[1]
    ni = b[2]
    nj = b[3]
    nk = b[4]
    t2 = N.zeros( (6,ni*nj*nk) )
    l = (0, 1, 2, 3, 4, 8)
    c = 0
    for i in l:
        t2[c] = t[i]
        c += 1

    v = ['x,y,z,ro,rovx,cellnf', t2, ni, nj, nk]
    out.append(v)

C.convertArrays2File(out, 'test.plt')

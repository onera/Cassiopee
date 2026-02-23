# - distribute (array) -
import Generator as G
import Distributor2 as D2
import numpy
from Converter.Internal import E_NpyInt

# Distribution sans communication entre blocs
N = 11
arrays = []
for i in range(N):
    a = G.cart( (0,0,0), (1,1,1), (10+i, 10, 10) )
    arrays.append(a)
out = D2.distribute(arrays, NProc=5); print(out)

# Distribution avec des perfos differentes pour chaque proc
out = D2.distribute(arrays, NProc=3, perfo=[(1,0,0), (1.2,0,0), (0.2,0,0)]); print(out)

# Distribution avec forcage du bloc 0 sur le proc 1, du bloc 2 sur le proc 3
# -1 signifie que le bloc est a equilibrer
prescribed = [-1 for x in range(N)]
prescribed[0] = 1; prescribed[2] = 3
out = D2.distribute(arrays, NProc=5, prescribed=prescribed); print(out)

# Distribution avec communications entre blocs, perfos identique pour tous
# les procs
volCom = numpy.zeros((N, N), dtype=E_NpyInt)
volCom[0,1] = 100; # Le bloc 0 echange 100 pts avec le bloc 1
out = D2.distribute(arrays, NProc=5, com=volCom, perfo=(1,0.,0.1)); print(out)

# Distribution avec des solveurs differents pour les blocs (le solveur est 2
# fois plus couteux pour les bloc 2 et 4)
out = D2.distribute(arrays, weight=[1,2,1,2,1,1,1,1,1,1,1], NProc=3); print(out)

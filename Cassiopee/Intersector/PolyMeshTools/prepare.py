import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import os, sys
import Intersector.PyTree as XOR

import os, sys
import Intersector.PyTree as XOR


if len(sys.argv) is not 3:
    print("ARG ERROR : 2 arguments to provide : mesh1 mesh2")
    sys.exit()

ifile1=sys.argv[1]
ifile2=sys.argv[2]

t1 = C.convertFile2PyTree(ifile1)
t2 = C.convertFile2PyTree(ifile2)

#print('getOverlappingFaces')
res = XOR.getOverlappingFaces(t1, t2, RTOL = 0.3, ps_min = 0.95)

# get pgids for t1 zones only : first par of each pairs
nb_zones = len(res)
t1zones_pgids = []
for i in range(nb_zones):
    t1zones_pgids.append(res[i][0])


#print('agglomerateCellsWithSpecifiedFaces')
XOR._agglomerateCellsWithSpecifiedFaces(t1, t1zones_pgids)

#print(t)
C.convertPyTree2File(t1, "out.cgns")

# - convertFile2Arrays (arrays) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

LOCAL = test.getLocal()

# Create test meshes
cart1 = G.cart((0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
cart2 = G.cartTetra((0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
cart3 = G.cartTetra((0,0,0), (0.1, 0.2, 1.), (11, 11, 1))
cart4 = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
l1 = D.line((0,0,50), (1,1,50))
l2 = D.line((0,0,1), (1,1,1))

# bin_tp
C.convertArrays2File([cart1, cart2], LOCAL+'/out.plt', 'bin_tp')
A = C.convertFile2Arrays(LOCAL+'/out.plt', 'bin_tp')
test.testA(A, 1)

# fmt_tp
C.convertArrays2File([cart1], LOCAL+'/out.tp', 'fmt_tp')
A = C.convertFile2Arrays(LOCAL+'/out.tp', 'fmt_tp')
test.testA(A, 2)

# fmt_tp
C.convertArrays2File([cart1,cart2,cart4], LOCAL+'/out.tp', 'fmt_tp')
A = C.convertFile2Arrays(LOCAL+'/out.tp', 'fmt_tp')
test.testA(A, 21)

# fmt_tp
C.convertArrays2File([cart1], LOCAL+'/out.tp', 'fmt_tp',zoneNames=['abcdefghijklmnopqrstuvwxyz'*200])
A = C.convertFile2Arrays(LOCAL+'/out.tp', 'fmt_tp')
test.testA(A, 22)

# bin_v3d
C.convertArrays2File([cart1], LOCAL+'/out.v3d', 'bin_v3d')
A = C.convertFile2Arrays(LOCAL+'/out.v3d', 'bin_v3d')
test.testA(A, 3)

# fmt_v3d
#C.convertArrays2File([cart1], LOCAL+'/out.dat', 'fmt_v3d')
#A = C.convertFile2Arrays(LOCAL+'/out.dat', 'fmt_v3d')
#test.testA(A, 4)

# bin_plot3d
#C.convertArrays2File([cart1], LOCAL+'/out.dat', 'bin_plot3d')
#A = C.convertFile2Arrays(LOCAL+'/out.dat', 'bin_plot3d')
#test.testA(A, 5)

# fmt_pov
C.convertArrays2File([cart3], LOCAL+'/out.pov', 'fmt_pov')
A = C.convertFile2Arrays(LOCAL+'/out.pov', 'fmt_pov')
test.testA(A, 6)

# fmt_mesh
C.convertArrays2File([cart3], LOCAL+'/out.mesh', 'fmt_mesh')
A = C.convertFile2Arrays(LOCAL+'/out.mesh', 'fmt_mesh')
test.testA(A, 7)

# fmt_su2
C.convertArrays2File([cart3], LOCAL+'/out.su2', 'fmt_su2')
A = C.convertFile2Arrays(LOCAL+'/out.su2', 'fmt_su2')
test.testA(A, 71)

# bin_stl
C.convertArrays2File([cart3], LOCAL+'/out.stl', 'bin_stl')
A = C.convertFile2Arrays(LOCAL+'/out.stl', 'bin_stl')
test.testA(A, 75)

# fmt_obj
C.convertArrays2File([cart3], LOCAL+'/out.obj', 'fmt_obj')
A = C.convertFile2Arrays(LOCAL+'/out.obj', 'fmt_obj')
test.testA(A, 77)

# bin_pickle
C.convertArrays2File([cart3], LOCAL+'/out.pickle', 'bin_pickle')
A = C.convertFile2Arrays(LOCAL+'/out.pickle', 'bin_pickle')
test.testA(A, 8)

# fmt_xfig
C.convertArrays2File([l1], LOCAL+'/out.fig', 'fmt_xfig')
A = C.convertFile2Arrays(LOCAL+'/out.fig', 'fmt_xfig')
test.testA(A, 11)

# fmt_svg
C.convertArrays2File([l2], LOCAL+'/out.svg', 'fmt_svg')
A = C.convertFile2Arrays(LOCAL+'/out.svg', 'fmt_svg')
test.testA(A, 12)

# fmt_svg avec densite
#C.convertArrays2File([l2], LOCAL+'/out.svg', 'fmt_svg')
#A = C.convertFile2Arrays(LOCAL+'/out.svg', 'fmt_svg', density=2.)
#test.testA(A, 13)

# fmt cedre
#C.convertArrays2File([cart4], LOCAL+'/out.d', 'fmt_cedre')
#A = C.convertFile2Arrays(LOCAL+'/out.d', 'fmt_cedre')
#test.testA(A, 14)

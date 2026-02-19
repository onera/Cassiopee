# - generateIBMMesh (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Post.PyTree as P
import Transform.PyTree as T
import Generator.IBM as G_IBM
import Geom.IBM as D_IBM
import KCore.test as test

#-------------------------------
# 2D
#-------------------------------
tb = D.circle((0.,0.,0.), 1., N=20)
tb = D_IBM.setDfar(tb, 10.)
tb = D_IBM.setSnear(tb, 0.1)
tb = C.newPyTree(['Base', tb])

xmin,xmax,ymin,ymax = 0.5, 9., -1, 1
tbox = D.box((xmin, ymin, 0), (xmax, ymax, 0), N=10, ntype='STRUCT')
tbox = P.exteriorFaces(tbox)
tbox = T.join(tbox)
tbox = D_IBM.setSnear(tbox, 0.1)
tbox = C.newPyTree(['Base', tbox])

t = G_IBM.generateIBMMesh(tb, dimPb=2, vmin=11, octreeMode=1, check=False)
test.testT(t,1)

t = G_IBM.generateIBMMesh(tb, dimPb=2, vmin=11, octreeMode=1, check=False, tbox=tbox)
test.testT(t,2)

#-------------------------------
# 3D
#-------------------------------
tb = D.sphere6((0.,0.,0.), 1., N=10, ntype='TRI')
tb = D_IBM.setDfar(tb, 10.)
tb = D_IBM.setSnear(tb, 0.1)
tb = C.newPyTree(['Base', tb])

xmin,xmax,ymin,ymax,zmin,zmax = 0.5, 9., -1, 1, -1, 1
tbox = D.box((xmin, ymin, zmin), (xmax, ymax, zmax), N=10, ntype='TRI')
tbox = D_IBM.setSnear(tbox, 0.1)
tbox = C.newPyTree(['Base', tbox])

t = G_IBM.generateIBMMesh(tb, dimPb=3, vmin=21, octreeMode=1, check=False)
test.testT(t,3)

t = G_IBM.generateIBMMesh(tb, dimPb=3, vmin=21, octreeMode=1, check=False, tbox=tbox)
test.testT(t,4)

#-------------------------------
# 3D + SYM
#-------------------------------
tb = D.cylinder((0,0,0), 0.25, 1, N=10, ntype='TRI')
tb = T.rotate(tb, (0,0,0), (1,0,0), 90)
tb = T.translate(tb, (0,1,0))
tb = D_IBM.setDfar(tb, 10.)
tb = D_IBM.setSnear(tb, 0.1)
tb = C.newPyTree(['Base', tb])
D_IBM._symetrizePb(tb, 'Base', snear_sym=1, dir_sym=2)

t = G_IBM.generateIBMMesh(tb, dimPb=3, vmin=21, octreeMode=1, check=False)
test.testT(t,5)
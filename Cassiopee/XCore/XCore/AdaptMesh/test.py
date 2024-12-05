import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import XCore.PyTree as X
import Converter.Internal as I
import Converter.elsAProfile as elsAProfile
import numpy as np
import Connector.PyTree as Conn
import Intersector.PyTree as XOR
import Post.PyTree as P

TOL_match = 1e-7
TOL_recover = 1e-11
rot_angle = 2.*np.pi/36.

# MASTER

m = C.convertFile2PyTree("data/maillage_veine.cgns")

(BCs_m, BCNames_m, BCTypes_m) = C.getBCs(m)
families_m = I.getNodesFromType2(m, 'Family_t')

m = C.convertArray2NGon(m)
m = T.merge(m)

# rotate master

rot_center = (0.,0.,0.)
rot_axis = (1.,0.,0.)

r1 = T.rotate(m, rot_center, rot_axis, rot_angle * 180. / np.pi)
r2 = T.rotate(m, rot_center, rot_axis, -rot_angle * 180. / np.pi)

# define helper tags

C._initVars(m, 'centers:keep', 1.0)
C._initVars(r1, 'centers:keep', 0.0)
C._initVars(r2, 'centers:keep', 0.0)

# merge

zm = I.getZones(m)[0]
zr1 = I.getZones(r1)[0]
zr2 = I.getZones(r2)[0]
mm = T.merge([zm, zr1, zr2])

# SLAVE

s = C.convertFile2PyTree("data/maillage_fente_axi.cgns")

(BCs_s, BCNames_s, BCTypes_s) = C.getBCs(s)
families_s = I.getNodesFromType2(s, 'Family_t')

# rotate slave

s1 = T.rotate(s, rot_center, rot_axis, -rot_angle * 180. / np.pi)

# REMOVE INTERSECTING PLANES

IM = X.IntersectMesh_Init(mm)
sp = X.removeIntersectingKPlanes(IM, s)
sp = T.merge(sp)
C._initVars(sp, 'centers:keep', 1.0)

sp1 = X.removeIntersectingKPlanes(IM, s1)
sp1 = T.merge(sp1)
C._initVars(sp1, 'centers:keep', 0.0)

zp = I.getZones(sp)[0]
zp1 = I.getZones(sp1)[0]
sp = T.merge([zp, zp1])
sfaces = X.extractFacesFromPointTag(sp, "tag")
print(len(sfaces))
cont = I.createUniqueChild(I.getZones(sp)[0], 'ZoneBC', 'ZoneBC_t')
I.newBC(name='CONTACT', pointList=sfaces[:]+1, btype='BCWall', parent=cont)
C.convertPyTree2File(sp, "sp.cgns")


# AdaptGeom

AM = X.AdaptMesh_Init(mm)

#sfaces = X.AdaptMesh_AdaptGeom(AM, sp, sfaces)
#cont = I.createUniqueChild(I.getZones(sp)[0], 'ZoneBC', 'ZoneBC_t')
#I.newBC(name='CONTACT', pointList=sfaces, btype='BCWall', parent=cont)
sa = X.AdaptMesh_AdaptGeom(AM, sp, sfaces)
C.convertPyTree2File(sa, "sa.cgns")

ma = X.AdaptMesh_ExtractMesh(AM)
C.convertPyTree2File(ma, "ma.cgns")

X.IntersectMesh_Exit(IM)
X.AdaptMesh_Exit(AM)

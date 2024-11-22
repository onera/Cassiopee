import OCC
import OCC.PyTree
import OCC.occ as occ
import Converter
import Converter.PyTree as C
import Generator
import Generator.PyTree as GP
import Converter.Internal as Internal
import Intersector.PyTree as XOR
import numpy as np
import sys
import KCore.test as test

test.TOLERANCE = 1.e-9

#import Ael.Quantum      as KDG

# get args
ifile='hammer.iges'
forma = 'fmt_iges'
height = 1
nlays = 1

hook = occ.readCAD(ifile, forma)
WALLBCS = ['BCWall', 'BCWallInviscid','BCWallViscous', 'BCWallViscousIsothermal', 'BCSymmetryPlane']


# DISCRETIZED
t = OCC.PyTree.convertCAD2PyTree(ifile, format=forma,
                          h=400, chordal_err=0, growth_ratio = 0., algo=1)


# EXTRUDED
# convert to NUGA format
t = XOR.convertTree2NUGANGON(t)

# set the surf as Wall BC
C._fillEmptyBCWith(t, 'wall', 'BCWall')
(BCs, BCNames, BCTypes) = C.getBCs(t)

# extrude
t = XOR.extrudeSurf(t, layer_height = height, nlayers = nlays, strategy = 1)

# put back BCs
C._recoverBCs(t, (BCs, BCNames, BCTypes), tol=1.e-6)


# Initialisation hx, hy, hz, ncadid
Internal.__FlowSolutionNodes__ = 'CADData'
t = C.initVars(t, 'hx', sys.float_info.max)
t = C.initVars(t, 'hy', sys.float_info.max)
t = C.initVars(t, 'hz', sys.float_info.max)
t = C.initVars(t, 'ncadid', -1)
# Initialisation u, v
t = C.initVars(t, 'u', sys.float_info.max)
t = C.initVars(t, 'v', sys.float_info.max)

# Initialisation fcadid
fcadid = np.empty((XOR.nb_faces(t)), dtype=np.int32)
fcadid[:] = -1

XOR._setZonesAndJoinsUId(t)
hmsh = XOR.createHMesh(t, 0)

#-------------------------------------------
# getNodalParameters
#-------------------------------------------
# Get zones
z = Internal.getZones(t)
# Get Information
c = C.getFields(Internal.__GridCoordinates__, z[0])
hx = C.getFields('CADData', C.extractVars(z, ['hx']))
hy = C.getFields('CADData', C.extractVars(z, ['hy']))
hz = C.getFields('CADData', C.extractVars(z, ['hz']))
ncadid =  C.getFields('CADData', C.extractVars(z, ['ncadid']))
u = C.getFields('CADData', C.extractVars(z, ['u']))
v = C.getFields('CADData', C.extractVars(z, ['v']))

# Get BC points
wall_face_ids = XOR.getBCPtListOfType(z, WALLBCS)

occ.getNodalParameters( c, wall_face_ids, hook, u, v, hx, hy, hz, ncadid)

# Mettre a jour des variables
t = C.setFields(hx,t, 'nodes')
t = C.setFields(hy,t, 'nodes')
t = C.setFields(hz,t, 'nodes')
t = C.setFields(ncadid,t, 'nodes')
t = C.setFields(u,t, 'nodes')
t = C.setFields(v,t, 'nodes')

#-------------------------------------------
# updateFcadidFromNcadid
#-------------------------------------------
p = Internal.getNodeFromName(t, 'CADData'); 
Internal.newDataArray('fcadid', value=fcadid, parent=p)

occ.updateFcadidFromNcadid( c, wall_face_ids, ncadid, fcadid)

#-------------------------------------------
# adaptCells
#-------------------------------------------
n = C.getNPts(t)
nv = np.empty((n,), dtype=np.int32)
nv[:] = 1

t = XOR.adaptCells(t,nv, sensor_type=2,hmesh = hmsh)
t = XOR.conformizeHMesh(t, hmsh)
t = XOR.closeCells(t)
Internal._rmNodesByName(t, 'zid')
Internal._rmNodesByName(t, 'rid')

#-------------------------------------------
# update updateNcadidFromFcadid
#-------------------------------------------
# Get input data
z = Internal.getZones(t)
c = C.getFields(Internal.__GridCoordinates__, z[0])
ncadid =  C.getFields('CADData', C.extractVars(z, ['ncadid']))
fcadid = Internal.getChildFromName(Internal.getNodeFromName(z, 'CADData'), 'fcadid')[1]
wall_face_ids = XOR.getBCPtListOfType(z, WALLBCS, None)

occ.updateNcadidFromFcadid(c, wall_face_ids, ncadid, fcadid)

# Update ncadid
t = C.setFields(ncadid, t, 'nodes')

#-------------------------------------------
# interpolateHMeshNodalField
#-------------------------------------------
t = XOR.interpolateHMeshNodalField(t, hmsh, ['u', 'v'])


#-------------------------------------------
# getNodalParameters
#-------------------------------------------
# Get zones
z = Internal.getZones(t)
# Get Information
c = C.getFields(Internal.__GridCoordinates__, z[0])
hx = C.getFields('CADData', C.extractVars(z, ['hx']))
hy = C.getFields('CADData', C.extractVars(z, ['hy']))
hz = C.getFields('CADData', C.extractVars(z, ['hz']))
ncadid =  C.getFields('CADData', C.extractVars(z, ['ncadid']))
u = C.getFields('CADData', C.extractVars(z, ['u']))
v = C.getFields('CADData', C.extractVars(z, ['v']))

occ.getNodalParameters(c, wall_face_ids, hook, u, v, hx, hy, hz, ncadid)


# Mettre a jour des variables
t = C.setFields(hx, t, 'nodes')
t = C.setFields(hy, t, 'nodes')
t = C.setFields(hz, t, 'nodes')
t = C.setFields(ncadid,t, 'nodes')
t = C.setFields(u, t, 'nodes')
t = C.setFields(v, t, 'nodes')

test.testT(t,1)
#C.convertPyTree2File(t, 'out.cgns')

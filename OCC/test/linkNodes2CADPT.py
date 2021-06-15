
import OCC
import OCC.PyTree
import OCC.occ as occ
import Converter
import Converter.PyTree as C
import Generator
import Generator.PyTree as G
import Converter.Internal as Internal
import Intersector.PyTree as XOR
import numpy, sys


# get args
ifile='hammer.iges'
forma = 'fmt_iges'
height = 0.1
nlays = 1
ofile='extruded.cgns'

hook = occ.readCAD(ifile, forma)
WALLBCS = ['BCWall', 'BCWallInviscid','BCWallViscous', 'BCWallViscousIsothermal', 'BCSymmetryPlane']


# DISCRETIZED
hmax = 0
chord_err = 0
t = OCC.PyTree.convertCAD2PyTree(ifile, format=forma,
                          h=hmax, chordal_err=chord_err, growth_ratio = 0., algo=1)
t = G.close(t, tol=2.5)
#C.convertPyTree2File(t, 'hamm.cgns')

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
#C.convertPyTree2File(t, 'layer.cgns')

#-------------------------------------------
# Link nodes2CAD
#-------------------------------------------

# Initialisation hx, hy, hz, ncadid
Internal.__FlowSolutionNodes__ = 'CADData'
t = C.initVars(t, 'hx', sys.float_info.max)
t = C.initVars(t, 'hy', sys.float_info.max)
t = C.initVars(t, 'hz', sys.float_info.max)
t = C.initVars(t, 'ncadid', -1)

# Initialisation fcadid
fcadid = numpy.empty((XOR.nb_faces(t)), dtype=numpy.int32)
fcadid[:] = -1

# adapt
zs = Internal.getZones(t)
n = C.getNPts(zs[0])
nv = numpy.empty((n,), dtype=numpy.int32)
nv[:] = 1

t = XOR.adaptCells(t,nv, sensor_type=2)
t = XOR.closeCells(t)
#C.convertPyTree2File(t, 'adapted.cgns')

# Get zones
z = Internal.getZones(t)
# Get Information
c = C.getFields(Internal.__GridCoordinates__, z[0])
hx = C.getFields('CADData', C.extractVars(z, ['hx']))
hy = C.getFields('CADData', C.extractVars(z, ['hy']))
hz = C.getFields('CADData', C.extractVars(z, ['hz']))
ncadid =  C.getFields('CADData', C.extractVars(z, ['ncadid']))

# Get BC points
wall_face_ids = XOR.getBCPtListOfType(z, WALLBCS)


occ.linkNodes2CAD( c, wall_face_ids, hook, hx, hy, hz, ncadid )

# Mettre a jour des variables
t = C.setFields(hx,t, 'nodes')
t = C.setFields(hy,t, 'nodes')
t = C.setFields(hz,t, 'nodes')
t = C.setFields(ncadid,t, 'nodes')

C.convertPyTree2File(t, 'out.cgns')

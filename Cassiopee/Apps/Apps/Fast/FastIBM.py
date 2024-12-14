from Generator.IBM import buildOctree, generateIBMMesh, createRefinementBodies

from Connector.IBM import prepareIBMData, dist2wallIBM, blankingIBM, buildFrontIBM, setInterpDataIBM, initializeIBM

from Geom.IBM import setSnear, _setSnear, setDfar, _setDfar, snearFactor, _snearFactor, setIBCType, changeIBCType, _changeIBCType, initOutflow, _initOutflow, initInj, _initInj, setFluidInside, setFluidOutside

from Generator.IBMmodelHeight import computeModelisationHeight, computeSnearOpt

from Apps.Fast.WindTunnelOutPres import getInfo, _setUpOutletPressure, getPointsFromTree, setupMachProbe, recordDataMach, _controlOutletPressureMachProbe

import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Converter.Filter as Filter
import Generator.PyTree as G
import Post.PyTree as P
import Geom.PyTree as D
import Post.Mpi as Pmpi
import Transform.PyTree as T
import Intersector.PyTree as XOR
import Post.Probe as Probe
import Distributor2.PyTree as D2

import numpy
import math
import os

# Ajouter l'initialisation de la pression, temperature et masse volumique
# à l'aide des formulations isentropiques (cf. tau de S Mouton)

###############
# Point Probes
###############

def createPointProbes(probe_in, probePointsList):
    listOfZones = []

    for cpt, (x,y,z) in enumerate(probePointsList):
        point = D.point((x,y,z))
        point[0] = "point_{:03d}".format(cpt)
        listOfZones.append(point)

    probe = C.newPyTree(['Base', listOfZones])

    if Cmpi.rank == 0: C.convertPyTree2File(probe, probe_in)

    Cmpi.barrier()

    return None

def initPointProbes(t, probe_in, fields, bufferSize=100, append=False, historyDirectory='.'):
    if isinstance(probe_in, str): probes = C.convertFile2PyTree(probe_in)
    else: probes = Internal.copyTree(probe_in)

    fields = ['centers:'+fname for fname in fields if 'centers' not in fname] #extraction from cell-centered t

    dictOfProbes = {}
    for z in Internal.getZones(probes):
        name = z[0]
        xnode = Internal.getNodeFromName(z,'CoordinateX')
        ynode = Internal.getNodeFromName(z,'CoordinateY')
        znode = Internal.getNodeFromName(z,'CoordinateZ')
        point = [Internal.getValue(xnode),Internal.getValue(ynode),Internal.getValue(znode)]
        dictOfProbes[name] = point

    if 'centers:Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    for pname in dictOfProbes.keys():
        point = dictOfProbes[pname]
        filename = "probe_{:s}.cgns".format(pname)
        filename = os.path.join(historyDirectory, filename)

        probe = Probe.Probe(filename, t, X=point, fields=fields, bufferSize=bufferSize, append=append)
        dictOfProbes[pname] = probe

    if 'centers:Mach' in fields: C._rmVars(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: C._rmVars(t, ['centers:Pressure'])

    Cmpi.barrier()

    return dictOfProbes

def _updatePointProbes(t, dictOfProbes, it, fields):
    fields = ['centers:'+fname for fname in fields if 'centers' not in fname] #extraction from cell-centered t

    if 'centers:Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    for name, probe in dictOfProbes.items():
        probe.extract(t, time=it)

    if 'centers:Mach' in fields: C._rmVars(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: C._rmVars(t, ['centers:Pressure'])

    Cmpi.barrier()

    return None

###############
# Surface Probes
###############

def generateIsoXSurface__(tb, x):
    bbox = G.bbox(tb)
    alpha= 0.05
    DY = bbox[4]-bbox[1]; DZ = bbox[5]-bbox[2]
    YMIN = bbox[1]-alpha*DY ; ZMIN = bbox[2]-alpha*DY
    LY = DY + 2*alpha*DY ; LZ = DZ + 2*alpha*DZ
    NJ = 51; NK = 51
    a = G.cart((x, YMIN, ZMIN), (1, LY/(NJ-1), LZ/(NK-1)), (1, NJ, NK))
    a = C.convertArray2Tetra(a)
    zones = Internal.getZones(tb)+[a]
    z = T.join(zones)
    z = XOR.conformUnstr(z, tol=1e-10, itermax=1)
    zones = T.splitManifold(z)
    candidates = []
    eps = 1e-4
    for z in zones:
        bboxz = G.bbox(z)
        if abs(bboxz[0]-x)<eps and abs(bboxz[3]-x)<eps:
            if bboxz[1]>=bbox[1]-eps and bboxz[4]<=bbox[4]+eps and bboxz[2]>=bbox[2]-eps and bboxz[5]<=bbox[5]+eps:
                candidates.append(z)
    for i,z in enumerate(candidates):
        z[0] = "zone_{:02d}".format(i)
    candidates = T.join(candidates)
    surface = C.newPyTree(["X_{:5.3f}".format(x), candidates])
    return surface

def createSurfaceProbes(tb, surface_in, probeSurfaceList):
    surfaces = []

    for x in probeSurfaceList:
        ts = generateIsoXSurface__(tb, x)
        G._getSmoothNormalMap(ts)
        surfaces.append(ts)

    ts = Internal.merge(surfaces)

    if Cmpi.size > 1:
        for b in Internal.getBases(ts):
            T._splitNParts(b, Cmpi.size)
            D2._distribute(b, Cmpi.size)

    if Cmpi.rank == 0: C.convertPyTree2File(ts, surface_in)

    Cmpi.barrier()

    return None

def initSurfaceProbes(t, tc, surface_in, fields, bufferSize=100, historyDirectory='.'):
    if isinstance(surface_in, str):
        if Cmpi.size > 1: probes = Cmpi.convertFile2PyTree(surface_in, proc=Cmpi.rank)
        else: probes = C.convertFile2PyTree(surface_in)
    else: probes = Internal.copyTree(surface_in)

    dictOfProbes = {}

    if 'Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    tcs = Internal.rmNodesFromType(tc, 'ZoneSubRegion_t')
    if Cmpi.size <= 1: Cmpi._setProc(tcs, 0) # Security for Probe functions
    for var in fields+['cellN']: C._cpVars(t, 'centers:'+var, tcs, var)

    for b in Internal.getBases(probes):
        sname = b[0]
        tbs_loc = C.newPyTree([sname, Internal.getZones(b)])

        filename = "surface_{:s}.cgns".format(sname)
        filename = os.path.join(historyDirectory, filename)

        probe = Probe.Probe(filename, tPermeable=tbs_loc, fields=fields, bufferSize=bufferSize)
        tcs_loc = probe.prepare(tcs)

        dictOfProbes[sname] = [probe, tbs_loc, tcs_loc]

    Cmpi.barrier()

    return dictOfProbes

def _updateSurfaceProbes(t, dictOfProbes, fields):
    if 'Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    for key in dictOfProbes:
        probe, tbs_loc, tcs_loc = dictOfProbes[key]
        for var in fields:
            C._cpVars(t, 'centers:'+var, tcs_loc, var)
            C._initVars(tbs_loc, var, 1)

        probe.extract(tcs_loc, time=1, onlyTransfer=True)

    if 'Mach' in fields: C._rmVars(t, ['centers:Mach'])
    if 'Pressure' in fields: C._rmVars(t, ['centers:Pressure'])

    Cmpi.barrier()

    return None

def getMassflow(t):
    C._initVars(t, '{massflow}={Density}*({sx}*{VelocityX}+{sy}*{VelocityY}+{sz}*{VelocityZ})')
    massflow = abs(Pmpi.integ(t, 'massflow')[0])
    return massflow

def integrateSurfaceProbes(dictOfProbes):
    massflows = []
    for key in dictOfProbes:
        probe, tbs_loc, tcs_loc = dictOfProbes[key]
        massflow_loc = getMassflow(tbs_loc) # intégration de la masse volumique sur chaque surface probe
        massflows.append(massflow_loc)

    Cmpi.barrier()

    return massflows

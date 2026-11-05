# run this in parallel with kpython
from FSDataManager import FSClac, FSMesh, FSError
import FSCGNSConverter.FSCGNSConverter as FSCGNSConverter

import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Converter.Mpi as Cmpi
import Connector.Mpi as Xmpi
import Ffx
import numpy
import time

fileName = 'case0.0001M.h5' # Non OK, Wall badly reconstructed
fileName = 'case0.001M.h5' # OK 2 cores
#fileName = 'case0.5M.h5' # OK 2 cores
#fileName = 'case1M.h5' # OK 2 cores
#fileName = 'case10M.h5' # OK 10 cores
#fileName = 'case30M.h5' # OK 30 cores
#fileName = 'case50M.h5' # OK 50 cores
#fileName = 'case100M.h5' # OK 96 cores

#=================================
# load and split with FSDM
#=================================
start0 = time.time()
clac = FSClac()
fsmesh = FSMesh(clac)

# global numbering
#fsmesh.ImportMeshHDF5(Filename='case.h5')

# local numbering
if Cmpi.master: print("Reading and splitting ...", flush=True)
start = time.time()
meshOps = FSCGNSConverter.buildMeshOps(fileName, partitioningLibrary="PARMETIS",
                                       preserveCellStacks=False,
                                       verbose=True)
fsmesh.DoOps(meshOps) or FSError.PrintAndExit()
if Cmpi.master: print("elapsed loadAndSplit: ", time.time()-start, flush=True)
start = time.time()

#=======================================
# convert to cgns with FSCGNSConverter
#=======================================
if Cmpi.master: print("Converting to CGNS ...", flush=True)
convObj = FSCGNSConverter.FSCGNSConverter(clac=clac,
                                          fsmesh=fsmesh,
                                          bcDict={1:'BCWall'},
                                          coordsName='Coordinates')
convObj.convert(forFFDX=True)
t = convObj.pyTree
#Cmpi.convertPyTree2File(t, 'out.cgns')
#exit()
if Cmpi.master: print("elapsed convert: ", time.time()-start, flush=True)
start = time.time()

# rename variables name
Internal._renameNode(t, "FlowSolution#CellAvgState", "FlowSolution#Centers")
#P._renameVars(t, ['centers:CellAvgState.Density'], ['centers:Density'])
#P._renameVars(t, ['centers:CellAvgState.MomentumX'], ['centers:MomentumX'])
#P._renameVars(t, ['centers:CellAvgState.MomentumY'], ['centers:MomentumY'])
#P._renameVars(t, ['centers:CellAvgState.MomentumZ'], ['centers:MomentumZ'])
#P._renameVars(t, ['centers:CellAvgState.EnergyStagnationDen'], ['centers:EnergyStagnationDensity'])
#P._renameVars(t, ['centers:CellAvgState.Viscosity_EddyMolec'], ['centers:Viscoisty_EddyMolecular'])
P._renameVars(t, ['centers:EnergyStagnationDen'], ['centers:EnergyStagnationDensity'])
P._renameVars(t, ['centers:Viscosity_EddyMolec'], ['centers:Viscosity_EddyMolecularRatio'])
#Internal._rmNodesFromName(t, 'FlowSolution#Centers')

#======================================================
# need to close for coherency between nodes and centers
#======================================================
#G._close(t)

#==============================
# convert2NGon avec recover BCs
#==============================
#if Cmpi.master: print("Converting to ngon ...", flush=True)
#C._convertArray2NGon(t, recoverBC=True)
#if Cmpi.master: print("elapsed convert2NGon: ", time.time()-start, flush=True)
#start = time.time()

#==============
# connectMatch
#==============
#if Cmpi.master: print("Connecting ...", flush=True)
zones = Internal.getZones(t)
Xmpi._connectMatchNGon(zones[0]) # connectMatchNGon limited to one zone per processor
if Cmpi.master: print("elapsed connectMatch: ", time.time()-start, flush=True)
start = time.time()

#=============
# Appel a FFX
#=============
state = {
    'Alpha': 2.0,
    'Roinf': 2.0,
    'Uinf': 1.0,
    'MomentumX': 2.0*numpy.cos(numpy.radians(2.)),
    'MomentumY': 0.,
    'MomentumZ': 2.0*numpy.sin(numpy.radians(2.)),
    'Lref': 50.,
    'Gamma': 1.4,
    'Pinf': 1.0,
    'Rgp': 0.03,
    'Csth': 110.4,
    'Tsth': 1.,
    'Musth': 1.716e-5,
    'Sref': 1.,
}

param = {
    'silent': True,
    'verbose': False
}

Cmpi.convertPyTree2File(t, 'out2.cgns')
exit()

if Cmpi.master: print("FFX ...", flush=True)
t, res = Ffx.computeOneAnalysis(t, state, param, xTP=2.0, filename=None)
if Cmpi.master: print("elapsed FFX: ", time.time()-start, flush=True)
start = time.time()

#===========
# sauvegarde
#===========
if Cmpi.master: print("Saving file ...", flush=True)


if Cmpi.master: print("full time elapsed: ", time.time()-start0, flush=True)

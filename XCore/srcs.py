import KCore.Dist as Dist
from KCore.config import *

SCOTCH=True; ZOLTAN=True
# 0: None, 1: paradigma, 2: paradigma23
PARADIGMA=2

(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                     additionalIncludePaths)

#==============================================================================
# Fichiers c++
#==============================================================================
cpp_srcs = ['XCore/CMP/src/recv_buffer.cpp', 
            'XCore/CMP/src/send_buffer.cpp',
            'XCore/xmpi/context_mpi_impl.cpp',
            'XCore/xmpi/context_stub_impl.cpp',
            'XCore/xmpi/communicator.cpp',
            'XCore/test/xmpi_t1.cpp',
            ]
if mpi: # source that requires mpi
    cpp_srcs += [
            'XCore/SplitElement/splitter.cpp',
            'XCore/chunk2partNGon.cpp',
            'XCore/chunk2partElt.cpp',
            'XCore/adaptMesh/cut.cpp',
            'XCore/adaptMesh/tree.cpp',
            'XCore/adaptMesh/mesh.cpp',
            'XCore/adaptMesh/math.cpp',
            'XCore/adaptMesh/adaptMesh.cpp',
            'XCore/adaptMesh/comm.cpp',
            'XCore/adaptMesh/metric.cpp',
            'XCore/adaptMesh/topo.cpp',
            'XCore/adaptMesh/distribute.cpp',
            'XCore/adaptMesh/adaptMeshSeq.cpp',
            'XCore/adaptMesh/createAdaptMesh.cpp',
            'XCore/adaptMesh/extractLeafMesh.cpp',
            'XCore/adaptMesh/adaptMeshDir.cpp',
            'XCore/adaptMesh/computeHessianNGon.cpp',
            'XCore/adaptMesh/computeGradientNGon.cpp',
            'XCore/adaptMesh/computeCellCentersAndVolumes.cpp',
            'XCore/common/mem.cpp',
            'XCore/exchangeFields.cpp'
            ]
else:
    cpp_srcs += [
        'XCore/chunk2partNGon_stub.cpp',
        'XCore/chunk2partElt_stub.cpp',
        'XCore/adaptMesh/adaptMesh_stub.cpp',
        'XCore/SplitElement/splitter_stub.cpp']

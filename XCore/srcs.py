import KCore.Dist as Dist
from KCore.config import *

SCOTCH=True; ZOLTAN=False
# 0: None, 1: paradigma, 2: paradigma23
PARADIGMA=0

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

            'XCore/exchangeFields.cpp',

            'XCore/common/mem.cpp',
            'XCore/common/common.cpp',

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
            

            'XCore/adaptMesh2/PrepareMeshForAdaptation.cpp',

            'XCore/adaptMesh2/AdaptMesh.cpp',
            'XCore/adaptMesh2/Mesh.cpp',
            'XCore/adaptMesh2/Topo.cpp',
            'XCore/adaptMesh2/Tree.cpp',
            'XCore/adaptMesh2/Refine.cpp',
            'XCore/adaptMesh2/Unrefine.cpp',
            'XCore/adaptMesh2/Geom.cpp',
            
            'XCore/adaptMesh2/Comm.cpp',
            
            'XCore/adaptMesh2/RenumberMesh.cpp',

            'XCore/adaptMesh2/Gradient.cpp',
            'XCore/adaptMesh2/Hessian.cpp',
            'XCore/adaptMesh2/HessianToMetric.cpp',
            'XCore/adaptMesh2/MetricToRefData.cpp',
            
            'XCore/adaptMesh2/IO.cpp',

            'XCore/adaptMesh2/Quad.cpp',
            'XCore/adaptMesh2/Tri.cpp',

            'XCore/adaptMesh2/Penta.cpp',
            'XCore/adaptMesh2/Hexa.cpp',
            'XCore/adaptMesh2/Tetra.cpp',
            'XCore/adaptMesh2/Pyra.cpp',

            'XCore/intersectMesh/removeIntersectingKPlanes.cpp',
            'XCore/intersectMesh/intersectSurf.cpp',
            'XCore/intersectMesh/smesh.cpp',

            'XCore/intersectMesh/mesh.cpp',
            'XCore/intersectMesh/io.cpp',
            'XCore/intersectMesh/bbox.cpp',
            'XCore/intersectMesh/geom.cpp',
            'XCore/intersectMesh/topo.cpp',
            
            'XCore/intersectMesh/dcel.cpp',
            'XCore/intersectMesh/event.cpp',
            'XCore/intersectMesh/hedge.cpp',
            'XCore/intersectMesh/predicate.cpp',
            'XCore/intersectMesh/segment.cpp',
            'XCore/intersectMesh/status.cpp',
            'XCore/intersectMesh/sweep.cpp',
            'XCore/intersectMesh/vertex.cpp'
            ]
else:
    cpp_srcs += [
        'XCore/chunk2partNGon_stub.cpp',
        'XCore/chunk2partElt_stub.cpp',
        'XCore/adaptMesh/adaptMesh_stub.cpp',
        'XCore/exchangeFields_stub',
        'XCore/SplitElement/splitter_stub.cpp']

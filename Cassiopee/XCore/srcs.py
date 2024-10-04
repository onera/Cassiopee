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
            
            'XCore/common/mem.cpp',
            'XCore/common/common.cpp',

            'XCore/intersectMesh/BVH.cpp',
            'XCore/intersectMesh/AABB.cpp',

            'XCore/intersectMesh/extract.cpp',

            'XCore/intersectMesh/IntersectMesh_Init.cpp',
            'XCore/intersectMesh/IntersectMesh_TriangulateFaceSet.cpp',
            'XCore/intersectMesh/IntersectMesh_ExtractMesh.cpp',
            'XCore/intersectMesh/IntersectMesh_Exit.cpp',
            'XCore/intersectMesh/DDA.cpp',

            'XCore/intersectMesh/intersectMesh.cpp',
            'XCore/intersectMesh/removeIntersectingKPlanes.cpp',
            'XCore/intersectMesh/prepareMeshesForIntersection.cpp',
            'XCore/intersectMesh/cycle.cpp',
            'XCore/intersectMesh/dcel.cpp',
            'XCore/intersectMesh/event.cpp',
            'XCore/intersectMesh/hedge.cpp',
            'XCore/intersectMesh/face.cpp',
            'XCore/intersectMesh/vertex.cpp',
            'XCore/intersectMesh/event.cpp',
            'XCore/intersectMesh/status.cpp',
            'XCore/intersectMesh/snode.cpp',
            'XCore/intersectMesh/queue.cpp',
            'XCore/intersectMesh/segment.cpp',
            'XCore/intersectMesh/mesh.cpp',
            'XCore/intersectMesh/meshRefine.cpp',
            'XCore/intersectMesh/meshTopo.cpp',
            'XCore/intersectMesh/smesh.cpp',
            'XCore/intersectMesh/smeshRefine.cpp',
            'XCore/intersectMesh/io.cpp',
            'XCore/intersectMesh/primitives.cpp',
            'XCore/intersectMesh/vec3.cpp',
            'XCore/intersectMesh/triangle.cpp',
            'XCore/intersectMesh/point.cpp',
            'XCore/intersectMesh/ray.cpp',
            'XCore/intersectMesh/karray.cpp',
            'XCore/intersectMesh/meshExport.cpp',


            'XCore/AdaptMesh/AdaptMesh_Init.cpp',
            'XCore/AdaptMesh/AdaptMesh_ExtractMesh.cpp',
            'XCore/AdaptMesh/AdaptMesh_AssignRefData.cpp',
            'XCore/AdaptMesh/AdaptMesh_Adapt.cpp',
            'XCore/AdaptMesh/AdaptMesh_Exit.cpp',

            'XCore/AdaptMesh/AdaptMesh_ExtractData.cpp',
            
            'XCore/AdaptMesh/MeshInit.cpp',
            'XCore/AdaptMesh/MeshOrient.cpp', 
            'XCore/AdaptMesh/MeshClean.cpp',
            'XCore/AdaptMesh/MeshTopo.cpp',
            'XCore/AdaptMesh/MeshIO.cpp',
            'XCore/AdaptMesh/MeshSmooth.cpp',

            'XCore/AdaptMesh/MeshConnectivity.cpp', 
            'XCore/AdaptMesh/MeshConformize.cpp', 
            'XCore/AdaptMesh/MeshRefine.cpp',
            'XCore/AdaptMesh/MeshIso.cpp',
            'XCore/AdaptMesh/MeshDir.cpp',

            'XCore/AdaptMesh/H27.cpp',
            'XCore/AdaptMesh/H18.cpp',
            'XCore/AdaptMesh/Tetra.cpp',
            'XCore/AdaptMesh/Penta.cpp',
            'XCore/AdaptMesh/Pyra.cpp',
            
            'XCore/AdaptMesh/Q9.cpp',
            'XCore/AdaptMesh/Q6.cpp',
            
            'XCore/AdaptMesh/Tri.cpp',

            'XCore/AdaptMesh/Edge.cpp',
            
            'XCore/AdaptMesh/Karray.cpp'
            ]
if mpi: # source that requires mpi
    cpp_srcs += [
            'XCore/SplitElement/splitter.cpp',
            'XCore/exchangeFields.cpp',
            'XCore/chunk2partNGon.cpp',
            'XCore/chunk2partElt.cpp',

            'XCore/AdaptMesh/AdaptMesh_LoadBalance.cpp',
            'XCore/AdaptMesh/MeshComm.cpp',
            ]
else:
    cpp_srcs += [
        'XCore/SplitElement/splitter_stub.cpp',
        'XCore/exchangeFields_stub.cpp',
        'XCore/chunk2partNGon_stub.cpp',
        'XCore/chunk2partElt_stub.cpp',
        'XCore/AdaptMesh/stubs/AdaptMesh_LoadBalance_stub.cpp'
        ]

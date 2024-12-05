# Test if libmpi exists ======================================================
import KCore.Dist as Dist
from KCore.config import *
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                     additionalIncludePaths)
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py(additionalLibPaths,
                                                        additionalIncludePaths)

#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ["Intersector/conformUnstr.cpp",

            "Intersector/booleanOperations.cpp",
            "Intersector/xcelln.cpp",
            "Intersector/selfX.cpp",
            "Intersector/P1ConservativeChimeraCoeffs.cpp",

            "PolyMeshTools/splitFaces.cpp",
            "PolyMeshTools/aggloFaces.cpp",
            "PolyMeshTools/aggloCells.cpp",
            "PolyMeshTools/splitCells.cpp",
            "PolyMeshTools/adaptCells.cpp",
            "Intersector/testm.cpp",
            "PolyMeshTools/utils.cpp"
            ]

if mpi and mpi4py:
    cpp_srcs += ["PolyMeshTools/adaptCells_mpi.cpp",
                 "PolyMeshTools/utils_mpi.cpp"]
else:
    cpp_srcs += ["PolyMeshTools/adaptCells_mpi_stub.cpp",
                 "PolyMeshTools/utils_mpi_stub.cpp"]

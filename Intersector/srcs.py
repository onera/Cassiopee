# Test if libmpi exists ======================================================
try:
    import KCore.Dist as Dist
    from KCore.config import *

    (mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                         additionalIncludePaths)
except:
    mpi = True

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
            "PolyMeshTools/utils.cpp",

            "Intersector/testm.cpp"
            
            ]

if mpi:
    cpp_srcs += ["PolyMeshTools/adaptCells_mpi.cpp"]
else:
    cpp_srcs += ["PolyMeshTools/adaptCells_mpi_stub.cpp"]


# =============================================================================
# THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
# =============================================================================
# File   : zipper/Make_tst.mk
# SVN    : $Rev$ $Date$
# Cksum  : 
# =============================================================================
#
E_PYINCLUDE1 = -I$(E_PPREFIX)/include/python$(E_PYVERSION) -I$(E_NUMPY)/numpy/core/include
E_PYINCLUDES = $(E_PYINCLUDE1) -I$(E_PPREFIX1)/include/python$(E_PYVERSION)
E_PYLIBS     = $(E_PPREFIX)/lib/python$(E_PYVERSION)/config
E_PYINC      = -I$(E_PYLIBS) $(E_PYINCLUDES)
E_PYLIB      = -L$(E_PYLIBS) -lpython$(E_PYVERSION)

# Object files to put into library
# 
E_LIBOBJLIST=\
Usurp.o\
usurp/usurp.o\
usurp/BuildGraphM.o usurp/CalculateForcesAndMomentsM.o\
usurp/CalculatePolygonAreaM.o \
usurp/CalculateProjectedAreasM.o usurp/CalculateVolumesM.o\
usurp/CheckInputM.o \
usurp/CheckPanelPairM.o usurp/CheckRatiosM.o usurp/CheckTrianglesM.o \
usurp/ConvertPanelToPolygonM.o usurp/CreateNewNodeListM.o\
usurp/CreateTrianglesM.o \
usurp/CutOutStringsM.o usurp/DeallocatePanelsM.o usurp/DefinePanelsM.o \
usurp/DetermineSurfaceOverlapM.o usurp/EdgeDataM.o\
usurp/ExaminePrioritiesM.o \
usurp/FoMoArraysM.o usurp/GarbageCollectorM.o usurp/GetUserInputM.o\
usurp/GraphProceduresM.o \
usurp/GroupInfoM.o usurp/InsertNodeM.o usurp/IntrTypeM.o\
usurp/LoadPanelsM.o \
usurp/LoadPanelsNPHASEM.o usurp/MapTriQM.o usurp/NetAllocM.o\
usurp/OutputForcesAndMomentsM.o \
usurp/OutputOVERFLOWM.o usurp/OverlappedM.o usurp/PatchInfoM.o\
usurp/PreReadBCinM.o \
usurp/PreReadBCSM.o usurp/PreReadGenericM.o usurp/PreReadNPHASEM.o\
usurp/PreReadOVERFLOWM.o \
usurp/PreReadUNCLEMM.o usurp/PrintPolyM.o usurp/PriPairsM.o\
usurp/ProcessPairProceduresM.o \
usurp/ReadPanelWeightsM.o usurp/ResolveEdgeIntersectionsM.o \
usurp/RotateBackM.o \
usurp/RTreeProceduresM.o usurp/RTreeTypesM.o usurp/ShrinkVertexListM.o\
usurp/sort2M.o \
usurp/StorePatchesCFDSHIPM.o usurp/StorePatchesGenericM.o\
usurp/StorePatchesUNCLEM.o \
usurp/StructuredPatchProceduresM.o\
usurp/TecplotNPHASEM.o usurp/TimeKeeperM.o \
usurp/TranslationArraysM.o usurp/TriQuadM.o usurp/TypesM.o\
usurp/UnusedVerticesM.o \
usurp/UserInputM.o usurp/VertexDataM.o usurp/VolumeArraysM.o\
usurp/WriteMixsurFMPM.o \
usurp/WritePanelWeightsM.o usurp/CommandLineM.o\
usurp/DPLRM.o usurp/my_cpu_timeM.o \
usurp/StorePatchesOVERFLOWM.o usurp/TecplotTrianglesM.o\
usurp/WriteGridIBIM.o \
usurp/WritePatchesM.o usurp/WriteTriQM.o usurp/ctype.o\
usurp/gpc.o usurp/gpc_c_translator_clip.o \
usurp/gpc_c_translator_strip.o usurp/share_eps.o usurp/triangle.o \
usurp/triangle_translator_c.o

OTHEROBJS=\
usurp/usurp.o\
usurp/BuildGraphM.o usurp/CalculateForcesAndMomentsM.o\
usurp/CalculatePolygonAreaM.o \
usurp/CalculateProjectedAreasM.o usurp/CalculateVolumesM.o\
usurp/CheckInputM.o \
usurp/CheckPanelPairM.o usurp/CheckRatiosM.o usurp/CheckTrianglesM.o \
usurp/ConvertPanelToPolygonM.o usurp/CreateNewNodeListM.o\
usurp/CreateTrianglesM.o \
usurp/CutOutStringsM.o usurp/DeallocatePanelsM.o usurp/DefinePanelsM.o \
usurp/DetermineSurfaceOverlapM.o usurp/EdgeDataM.o\
usurp/ExaminePrioritiesM.o \
usurp/FoMoArraysM.o usurp/GarbageCollectorM.o usurp/GetUserInputM.o\
usurp/GraphProceduresM.o \
usurp/GroupInfoM.o usurp/InsertNodeM.o usurp/IntrTypeM.o\
usurp/LoadPanelsM.o \
usurp/LoadPanelsNPHASEM.o usurp/MapTriQM.o usurp/NetAllocM.o\
usurp/OutputForcesAndMomentsM.o \
usurp/OutputOVERFLOWM.o usurp/OverlappedM.o usurp/PatchInfoM.o\
usurp/PreReadBCinM.o \
usurp/PreReadBCSM.o usurp/PreReadGenericM.o usurp/PreReadNPHASEM.o\
usurp/PreReadOVERFLOWM.o \
usurp/PreReadUNCLEMM.o usurp/PrintPolyM.o usurp/PriPairsM.o\
usurp/ProcessPairProceduresM.o \
usurp/ReadPanelWeightsM.o usurp/ResolveEdgeIntersectionsM.o \
usurp/RotateBackM.o \
usurp/RTreeProceduresM.o usurp/RTreeTypesM.o usurp/ShrinkVertexListM.o\
usurp/sort2M.o \
usurp/StorePatchesCFDSHIPM.o usurp/StorePatchesGenericM.o\
usurp/StorePatchesUNCLEM.o \
usurp/StructuredPatchProceduresM.o\
usurp/TecplotNPHASEM.o usurp/TimeKeeperM.o \
usurp/TranslationArraysM.o usurp/TriQuadM.o usurp/TypesM.o\
usurp/UnusedVerticesM.o \
usurp/UserInputM.o usurp/VertexDataM.o usurp/VolumeArraysM.o\
usurp/WriteMixsurFMPM.o \
usurp/WritePanelWeightsM.o usurp/CommandLineM.o\
usurp/DPLRM.o usurp/my_cpu_timeM.o \
usurp/StorePatchesOVERFLOWM.o usurp/TecplotTrianglesM.o\
usurp/WriteGridIBIM.o \
usurp/WritePatchesM.o usurp/WriteTriQM.o usurp/ctype.o\
usurp/gpc.o usurp/gpc_c_translator_clip.o \
usurp/gpc_c_translator_strip.o usurp/share_eps.o usurp/triangle.o \
usurp/triangle_translator_c.o

E_LIB=\
Usurpl.so
#
# --------------------------------------------------------LAST LINE------

#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ['Connector/KInterp/BlkInterp.cpp',
            'Connector/KInterp/KMesh.cpp',
            'Connector/KInterp/BlkIntTreeNode.cpp',
            'Connector/KInterp/BlkInterpData.cpp',
            'Connector/KInterp/BlkInterpDataStruct.cpp',
            'Connector/KInterp/BlkInterpWithKMesh.cpp',
            'Connector/KInterp/BlkInterpAdt.cpp',
            'Connector/KInterp/BlkInterpAdt_getCell.cpp',
            'Connector/optimizeOverlap.cpp',
            'Connector/maximizeBlankedCells.cpp',
            'Connector/maskXRay.cpp',
            'Connector/maskGen.cpp',
            'Connector/blankCells.cpp',
            'Connector/blankCellsTetra.cpp',
            'Connector/getIntersectingDomainsAABB.cpp',
            'Connector/changeWall.cpp',
            'Connector/changeWallEX.cpp',
            'Connector/setDoublyDefinedBC.cpp',
            'Connector/getInterpolatedPoints.cpp',
            'Connector/getInterpolatedEXPoints.cpp',
            'Connector/setInterpolations.cpp',
            'Connector/setInterpData.cpp',
            'Connector/setInterpDataCons.cpp',
            'Connector/setInterpDataForGhostCells.cpp',
            'Connector/setInterpDataForGhostCellsNGon.cpp',
            'Connector/initNuma.cpp',
            'Connector/setInterpTransfers.cpp',
            'Connector/setInterpTransfersD.cpp',
            'Connector/IBC/setIBCTransfersD.cpp',
            'Connector/IBC/blankClosestTargetCells.cpp',
            'Connector/writeCoefs.cpp',
            'Connector/chimeraTransfer.cpp',
            'Connector/transferVariables.cpp',
            'Connector/blankIntersectingCells.cpp',
            'Connector/cellN2OversetHoles.cpp',
            'Connector/identifyMatching.cpp',
            'Connector/identifyDegenerated.cpp',
            'Connector/gatherMatching.cpp',
            'Connector/gatherMatchingNM.cpp',
            'Connector/gatherMatchingNGon.cpp',
            'Connector/gatherDegenerated.cpp',
            'Connector/IBC/setIBCTransfers.cpp',
            'Connector/IBC/setInterpDataMLS_IBMWall.cpp',
            'Connector/setInterpDataLS.cpp',
            'Connector/modifyBorders.cpp',
            'Connector/applyBCOverlaps.cpp',
            "Connector/getExtrapAbsCoefs.cpp",
            "Connector/getEmptyBCInfoNGON.cpp",
            "Connector/updateNatureForIBM.cpp",
            "Connector/getIBMPtsWithFront.cpp",
            "Connector/getIBMPtsWithTwoFronts.cpp",
            "Connector/getIBMPtsWithoutFront.cpp",
            "Connector/getIBMPtsBasic.cpp",
            "Connector/indiceToCoord2.cpp",
            "Connector/correctCoeffList.cpp",
            "Connector/modCellN.cpp",
            "Connector/LBM/setInterpTransfers.cpp",
            "Connector/IBC/LBM/setIBCTransfers.cpp",
            "Connector/computeFrictionVelocityIBM.cpp"]

#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ['Connector/Fortran/spalart_1d.for', # called in setIBCTRansfers.cpp/connector.h (spalart_1d_)
            'Connector/Fortran/BlkAdjustCellNatureFieldF.for', # called in blankCells.cpp (k6adjustcellnaturefield_)
            'Connector/Fortran/MaskSearchBlankedNodesXF.for', # called in blankCells.cpp (k6searchblankednodesx_)
            'Connector/Fortran/MaskSearchBlankedNodesX2DF.for', # called in blankCells.cpp (k6searchblankednodesx2d)
            'Connector/Fortran/MaskSearchBlankedNodesXDF.for', # called in blankCells.cpp (k6searchblankednodesxd)
            'Connector/Fortran/MaskSearchBlankedNodesXD2DF.for', # called in blankCells.cpp (k6searchblankednodesxd2d_)
            'Connector/Fortran/MaskSearchBlankedCellsX2DF.for', # called in blankCells.cpp (k6searchblankedcellsx2d_)
            'Connector/Fortran/MaskSearchBlankedCellsXF.for', # called in blankCells.cpp (k6searchblankedcellsx_)
            'Connector/Fortran/MaskSearchBlankedCellsX12DF.for', # called in blankCells.cpp (k6searchblankedcellsx12d_)
            'Connector/Fortran/MaskSearchBlankedCellsX1F.for', # called in blankCells.cpp (k6searchblankedcellsx1_)
            'Connector/Fortran/MaskSearchBlankedCellsX22DF.for', # called in blankCells.cpp (k6searchblankedcellsx22d_)
            'Connector/Fortran/MaskSearchBlankedCellsX2F.for', # called in blankCells.cpp (k6searchblankedcellsx2_)
            'Connector/Fortran/MaskSearchBlankedCellsXDF.for', # called in blankCells.cpp (k6searchblankedcellsxd_)
            'Connector/Fortran/MaskSearchBlankedCellsXD2DF.for', # called in blankCells.cpp (k6searchblankedcellsxd2d_)
            'Connector/Fortran/MaskSearchBlankedCellsXD2F.for', # called in blankCells.cpp (k6searchblankedcellsxd2_)
            'Connector/Fortran/MaskSearchBlankedCellsXD22DF.for', # called in blankCells.cpp (k6searchblankedcellsxd22d_)
            'Connector/Fortran/MaskSearchBlankedCellsXD1F.for', # called in blankCells.cpp (k6searchblankedcellsxd1_)
            'Connector/Fortran/MaskSearchBlankedCellsXD12DF.for', # called in blankCells.cpp (k6searchblankedcellsxd12d_)
            'Connector/Fortran/MaskSearchBlankedCellsUnstrXF.for', # called in blankCells.cpp (k6searchblankedcellstrix_)
            'Connector/Fortran/MaskSearchBlankedCellsUnstrXDF.for', # called in blankCells.cpp (k6searchblankedcellstrixd_)
            'Connector/Fortran/MaskSearchBlankedCellsUnstrX2F.for'] # called in blankCells.cpp (k6searchblankedcellstrix2_)

"""Toolbox for IBM preprocessing"""

# This module is deprecated, please use Connector.IBM instead

from Connector.IBM import _blankClosestTargetCells, _removeBlankedGrids, blankByIBCBodies, getIBMFront, getIBMFrontType1__, getIBMFrontType0__, _pushBackImageFront2, \
    _smoothImageFront, gatherFront, doInterp, _extractIBMInfo_param, \
    extractIBMInfo, getAllIBMPoints, prepareIBMData_legacy, createWallAdapt, \
    createIBMWZones, _computeKcurvParameter, _signDistance

from Generator.IBM import generateCartMesh__, adaptIBMMesh, generateIBMMesh, buildOctree, addRefinementZones__, octree2StructLoc__, \
    mergeByParent__, buildParentOctrees__, _addBCOverlaps, _addExternalBCs, \
    _modifPhysicalBCs__

from Generator.IBMmodelHeight import getMinimumCartesianSpacing, computeYplus, computeModelisationHeight, computeBestModelisationHeight

from Post.IBM import extractIBMWallFields

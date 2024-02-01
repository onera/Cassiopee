"""Toolbox for IBM preprocessing"""

# This module is deprecated, please use Connector.IBM instead

from Connector.IBM import _blankClosestTargetCells, _removeBlankedGrids, blankByIBCBodies, _addBCOverlaps, _addExternalBCs, \
    _modifPhysicalBCs__, getIBMFront, getIBMFrontType1__, getIBMFrontType0__, getIBMFrontType0Old__, _pushBackImageFront2, \
    _smoothImageFront, gatherFront, doInterp, doInterp2, doInterp3, _extractIBMInfo_param, \
    extractIBMInfo, getAllIBMPoints, prepareIBMData, prepareIBMData2, createWallAdapt, \
    createIBMWZones, _computeKcurvParameter, _signDistance

from Generator.IBM import generateCartMesh__, adaptIBMMesh, generateIBMMesh, buildOctree, addRefinementZones, octree2StructLoc__, \
    mergeByParent__, buildParentOctrees__

from Generator.IBMmodelHeight import getMinimumCartesianSpacing, computeYplus, computeModelisationHeight, computeBestModelisationHeight

from Post.IBM import extractIBMWallFields

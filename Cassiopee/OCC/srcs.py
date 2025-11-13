GORDON = True

import glob
import os

exceptions = ['../../OCC/occ_src/TKXSBase/stepread.c',
              '../../OCC/occ_src/TKXSBase/step.tab.c',
              '../../OCC/occ_src/TKXSBase/lex.step.c']

def getFiles(module):
    if os.path.exists('../../OCC/occ_src/%s/PACKAGES'%module):
        f = open('../../OCC/occ_src/%s/PACKAGES'%module, 'r')
        lines = f.readlines()
        f.close()
        dirs = []
        for l in lines:
            l = l[:-1]
            dirs.append('../../OCC/occ_src/%s'%l)
    else: dirs = ['../../OCC/occ_src/%s'%module]
    files = []
    for d in dirs:
        files += glob.glob(d+'/*.cxx')
        files += glob.glob(d+'/*.c')
    srcs = []
    for f in files:
        if f not in exceptions: srcs.append(f.replace('../../',''))
    return srcs

srcs = ['OCC/import_OCC_CAD_wrapper.cpp',
        'OCC/CADviaOCC.cpp',
        'OCC/OCCSurface.cpp',
        'OCC/CADviaOCC2.cpp',
        'OCC/OCCSurface2.cpp',
        'OCC/convertCAD2Arrays0.cpp',
        'OCC/convertCAD2Arrays1.cpp',
        'OCC/convertCAD2Arrays2.cpp',

        'OCC/Atomic/readCAD.cpp',
        'OCC/Atomic/writeCAD.cpp',
        'OCC/Atomic/createEmptyCAD.cpp',
        'OCC/Atomic/mergeCAD.cpp',
        'OCC/Atomic/freeHook.cpp',

        'OCC/Atomic/printOCAF.cpp',
        'OCC/Atomic/getFaceNameInOCAF.cpp',

        'OCC/Atomic/bottle.cpp',
        'OCC/Atomic/addSphere.cpp',
        'OCC/Atomic/addCylinder.cpp',
        'OCC/Atomic/addBox.cpp',
        'OCC/Atomic/addLine.cpp',
        'OCC/Atomic/addCircle.cpp',
        'OCC/Atomic/addSquare.cpp',
        'OCC/Atomic/addSpline.cpp',
        'OCC/Atomic/addArc.cpp',
        'OCC/Atomic/addGordonSurface.cpp',

        'OCC/Atomic/meshEdge.cpp',
        'OCC/Atomic/meshEdge2.cpp',
        'OCC/Atomic/identifyLoopsInEdges.cpp',
        'OCC/Atomic/parameterEdges.cpp',
        'OCC/Atomic/evalEdge.cpp',
        'OCC/Atomic/evalFace.cpp',
        'OCC/Atomic/projectOnFace.cpp',
        'OCC/Atomic/projectOnEdge.cpp',
        'OCC/Atomic/linkNodes2CAD.cpp',
        'OCC/Atomic/trimesh.cpp',
        'OCC/Atomic/analyse.cpp',
        'OCC/Atomic/getFaceArea.cpp',
        'OCC/Atomic/areEdgeIdentical.cpp',

        'OCC/Atomic/splitter.cpp',
        'OCC/Atomic/fix.cpp',
        'OCC/Atomic/trim.cpp',
        'OCC/Atomic/sewing.cpp',
        'OCC/Atomic/removeFaces.cpp',
        'OCC/Atomic/fillHole.cpp',
        'OCC/Atomic/addFillet.cpp',
        'OCC/Atomic/mergeFaces.cpp',
        'OCC/Atomic/loft.cpp',
        'OCC/Atomic/revolve.cpp',

        'OCC/Atomic/translate.cpp',
        'OCC/Atomic/scale.cpp',
        'OCC/Atomic/rotate.cpp',

        'OCC/Atomic/intersectEdgeFace.cpp',

        'OCC/Atomic/getOppData.cpp',
        'OCC/Atomic/identifyTags.cpp']

if GORDON:
    srcs += [
        'OCC/Gordon/BSplineAlgorithms.cpp',
        'OCC/Gordon/CurveNetworkSorter.cpp',
        'OCC/Gordon/Error.cpp',
        'OCC/Gordon/InterpolateCurveNetwork.cpp',
        'OCC/Gordon/PointsToBSplineInterpolation.cpp',
        'OCC/Gordon/BSplineApproxInterp.cpp',
        'OCC/Gordon/CurvesToSurface.cpp',
        'OCC/Gordon/GordonSurfaceBuilder.cpp',
        'OCC/Gordon/IntersectBSplines.cpp',
        'OCC/Gordon/occ_gordon.cpp']

#====================================================================================
import KCore.Dist as Dist
allMods = Dist.getOCCModules()

mod_srcs = {}
for m in allMods:
    mod_srcs[m] = getFiles(m)

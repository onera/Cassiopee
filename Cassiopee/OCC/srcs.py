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

        'OCC/Atomic/bottle.cpp',
        'OCC/Atomic/addSphere.cpp',
        'OCC/Atomic/addSquare.cpp',
        'OCC/Atomic/addCylinder.cpp',
        'OCC/Atomic/addBox.cpp',

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

        'OCC/Atomic/translate.cpp',
        'OCC/Atomic/scale.cpp',
        'OCC/Atomic/rotate.cpp',

        'OCC/Atomic/getOppData.cpp']

import KCore.Dist as Dist

allMods6 = ["FWOSPlugin", "TKPLCAF", "PTKernel", "TKPrim",
            "TKAdvTools", "TKPShape", "TKBinL", "TKService",
            "TKBin", "TKShapeSchema", "TKBinTObj", "TKShHealing",
            "TKBinXCAF", "TKStdLSchema", "TKBool", "TKStdSchema",
            "TKBO", "TKSTEP209", "TKBRep", "TKSTEPAttr",
            "TKCAF", "TKSTEPBase", "TKCDF", "TKSTEP",
            "TKernel", "TKSTL", "TKFeat", "TKTObj",
            "TKFillet", "TKTopAlgo", "TKG2d", "TKV3d",
            "TKG3d", "TKVoxel", "TKGeomAlgo", "TKVRML",
            "TKGeomBase", "TKXCAFSchema", "TKHLR", "TKXCAF",
            "TKIGES", "TKXDEIGES", "TKLCAF", "TKXDESTEP",
            "TKMath", "TKXMesh", "TKMesh", "TKXmlL",
            "TKMeshVS", "TKXml", "TKNIS", "TKXmlTObj",
            "TKOffset", "TKXmlXCAF", "TKOpenGl", "TKXSBase",
            "TKPCAF"]

allMods75 = ["TKBinL", "TKBin", "TKBinTObj",
             "TKBinXCAF", "TKBool", "TKBO", "TKBRep",
             "TKCAF", "TKCDF", "TKDCAF", "TKDraw", "TKernel",
             "TKFeat", "TKFillet", "TKG2d", "TKG3d", "TKGeomAlgo",
             "TKGeomBase", "TKHLR", "TKIGES", "TKLCAF", "TKMath",
             "TKMesh", "TKMeshVS", "TKOffset", "TKOpenGl", "TKPrim",
             "TKQADraw", "TKRWMesh", "TKService", "TKShHealing", "TKStdL",
             "TKStd", "TKSTEP209", "TKSTEPAttr", "TKSTEPBase", "TKSTEP",
             "TKSTL", "TKTObjDRAW", "TKTObj", "TKTopAlgo", "TKTopTest",
             "TKV3d", "TKVCAF", "TKViewerTest", "TKVRML", "TKXCAF", "TKXDEDRAW",
             "TKXDEIGES", "TKXDESTEP", "TKXMesh", "TKXmlL", "TKXml", "TKXmlTObj",
             "TKXmlXCAF", "TKXSBase", "TKXSDRAW"]

allMods75W = ["TKBin", "TKBinL", "TKBinTObj", "TKBinXCAF", "TKBO",
              "TKBool", "TKBRep", "TKCAF", "TKCDF", "TKernel",
              "TKFeat", "TKFillet", "TKG2d", "TKG3d", "TKGeomAlgo",
              "TKGeomBase", "TKHLR", "TKIGES", "TKLCAF", "TKMath",
              "TKMesh", "TKMeshVS", "TKOffset", "TKOpenGl",
              "TKPrim", "TKService",
              "TKShHealing", "TKStdLSchema",
              "TKStdSchema", "TKSTEP", "TKSTEP209", "TKSTEPAttr",
              "TKSTEPBase", "TKSTL", "TKTObj", "TKTopAlgo",
              "TKV3d", "TKVoxel", "TKVRML", "TKXCAF", "TKXCAFSchema",
              "TKXDEIGES", "TKXDESTEP", "TKXMesh", "TKXml",
              "TKXmlL", "TKXmlTObj", "TKXmlXCAF", "TKXSBase",
              "TKPCAF", "TKPLCAF", "TKNIS", "TKPShape", "TKShapeSchema"]


if Dist.getSystem()[0] == 'mingw': allMods = allMods75W
else: allMods = allMods75

mod_srcs = {}
for m in allMods:
    mod_srcs[m] = getFiles(m)

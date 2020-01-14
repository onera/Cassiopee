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

srcs = ['OCC/pyOCC.cpp', 
        'OCC/import_OCC_CAD_wrapper.cpp',
        'OCC/CADviaOCC.cpp', 
        'OCC/OCCSurface.cpp',
        'OCC/convertCAD2Arrays.cpp']

allMods = ["TKBin", "TKBinL", "TKBinTObj", "TKBinXCAF", "TKBO",
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

mod_srcs = {}
for m in allMods:
    mod_srcs[m] = getFiles(m)

#TKBin_srcs = getFiles("TKBin")
#TKBinL_srcs = getFiles("TKBinL")
#TKBinTObj_srcs = getFiles("TKBinTObj")
#TKBinXCAF_srcs = getFiles("TKBinXCAF")
#TKBO_srcs = getFiles("TKBO")
#TKBool_srcs = getFiles("TKBool")
#TKBrep_srcs = getFiles("TKBRep")
#TKCAF_srcs = getFiles("TKCAF")
#TKCDF_srcs = getFiles("TKCDF")
#TKernel_srcs = getFiles("TKernel"),
#TKFeat_srcs = getFiles("TKFeat")
#TKFillet_srcs = getFiles("TKFillet")
#TKG2d_srcs = getFiles("TKG2d")
#TKG3d_srcs = getFiles("TKG3d") 
#TKGeomAlgo_srcs = getFiles("TKGeomAlgo")
#TKGeomBase_srcs = getFiles("TKGeomBase")
#TKHLR_srcs = getFiles("TKHLR")
#TKIGES_srcs = getFiles("TKIGES")
#TKLCAF_srcs = getFiles("TKLCAF"), 
#TKMath_srcs = getFiles("TKMath"),
#TKMesh_srcs = getFiles("TKMesh")
#TKMeshVS_srcs = getFiles("TKMeshVS")
#TKNIS_srcs = getFiles("TKNIS")
#TKOffset_srcs = getFiles("TKOffset"), 
#TKOpenGl_srcs = getFiles("TKOpenGl") 
#TKPCAF_srcs = getFiles("TKPCAF")
#TKPLCAF_srcs = getFiles("TKPLCAF")
#TKPrim_srcs = getFiles("TKPrim")
#TKPShape_srcs = getFiles("TKPShape"), 
#TKService_srcs = getFiles("TKService"),
#TKShapeSchema_srcs = getFiles("TKShapeSchema")
#TKShHealing_srcs = getFiles("TKShHealing")
#TKStdLSchema_srcs = getFiles("TKStdLSchema"),
#TKStdSchema_srcs = getFiles("TKStdSchema")
#TKSTEP_srcs = getFiles("TKSTEP")
#TKSTEP209_srcs = getFiles("TKSTEP209")
#TKSTEPAttr_srcs = getFiles("TKSTEPAttr"),
#TKSTEPBase_srcs = getFiles("TKSTEPBase") 
#TKSTL_srcs = getFiles("TKSTL")
#TKTObj_srcs = getFiles("TKTObj")
#TKTopAlgo_srcs = getFiles("TKTopAlgo"),
#TKV3d_srcs = getFiles("TKV3d")
#TKVoxel_srcs = getFiles("TKVoxel")
#TKVRML_srcs = getFiles("TKVRML")
#TKXCAF_srcs = getFiles("TKXCAF"), 
#TKCAFSchema_srcs = getFiles("TKXCAFSchema"),
#TKDEIGES_srcs = getFiles("TKXDEIGES")
#TKXDESTEP_srcs = getFiles("TKXDESTEP")
#TKXMesh_srcs = getFiles("TKXMesh")
#TKXml_srcs = getFiles("TKXml"),
#TKXmlL_srcs = getFiles("TKXmlL")
#TKXmlTObj_srcs = getFiles("TKXmlTObj")
#TKXmlXCAF_srcs = getFiles("TKXmlXCAF")
#TKXSBase_srcs = getFiles("TKXSBase")

# Resticted lib
#TKernel_srcs = getFiles('TKernel')
#TKMath_srcs = getFiles('TKMath')
#TKGeomBase_srcs = getFiles('TKGeomBase')
#TKG2d_srcs = getFiles('TKG2d')
#TKG3d_srcs = getFiles('TKG3d')
#TKBrep_srcs = getFiles('TKBRep')
#TKGeomAlgo_srcs = getFiles('TKGeomAlgo')
#TKBool_srcs = getFiles('TKBool')
#TKPrim_srcs = getFiles('TKPrim')
#TKShHealing_srcs = getFiles('TKShHealing')
#TKTopAlgo_srcs = getFiles('TKTopAlgo')
#TKXSBase_srcs = getFiles('TKXSBase')
#TKIGES_srcs = getFiles('TKIGES')
#TKSTEP_srcs = getFiles('TKSTEP')
#TKSTEP2_srcs = getFiles('TKSTEP2')

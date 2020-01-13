import glob

exceptions = ['../../OCC/occ_src/TKXSBase/stepread.c', 
              '../../OCC/occ_src/TKXSBase/step.tab.c',
              '../../OCC/occ_src/TKXSBase/lex.step.c']

def getFiles(module):
    files = glob.glob('../../OCC/occ_src/%s/*.cxx'%module)
    files += glob.glob('../../OCC/occ_src/%s/*.c'%module)
    srcs = []
    for f in files:
        if f not in exceptions: srcs.append(f.replace('../../',''))
    return srcs

srcs = ['OCC/pyOCC.cpp', 
        'OCC/import_OCC_CAD_wrapper.cpp',
        'OCC/CADviaOCC.cpp', 
        'OCC/OCCSurface.cpp',
        'OCC/convertCAD2Arrays.cpp']

# TKernel
TKernel_srcs = getFiles('TKernel')

# TKMath
TKMath_srcs = getFiles('TKMath')

# TKGeomBase
TKGeomBase_srcs = getFiles('TKGeomBase')

# TKG2d
TKG2d_srcs = getFiles('TKG2d')

# TKG3d
TKG3d_srcs = getFiles('TKG3d')

# TKBRep
files = getFiles('TKBRep')
TKBRep_srcs = files[0:100]
TKBRep2_srcs = files[100:]

# TKGeomAlgo
TKGeomAlgo_srcs = getFiles('TKGeomAlgo')

# TKBool
files = getFiles('TKBool')
TKBool_srcs = files[0:100]
TKBool2_srcs = files[100:200]
TKBool3_srcs = files[200:300]
TKBool4_srcs = files[300:]

# TKPrim
TKPrim_srcs = getFiles('TKPrim')

# TKShHealing
TKShHealing_srcs = getFiles('TKShHealing')

# TKTopAlgo
TKTopAlgo_srcs = getFiles('TKTopAlgo')

# TKXSBase
TKXSBase_srcs = getFiles('TKXSBase')

# TKIGES
TKIGES_srcs = getFiles('TKIGES')

# TKSTEP
TKSTEP_srcs = getFiles('TKSTEP')

# TKSTEP2
TKSTEP2_srcs = getFiles('TKSTEP2')

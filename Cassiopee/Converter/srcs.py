EXPRESSION = True

try:
    import KCore.Dist as Dist
    from KCore.config import *
    (hdf, hdfIncDir, hdfLibDir, hdflibs) = Dist.checkHdf(additionalLibPaths, additionalIncludePaths)
    (netcdf, netcdfIncDir, netcdfLibDir, netcdflibs) = Dist.checkNetcdf(additionalLibPaths, additionalIncludePaths)

except ModuleNotFoundError:
    hdf = True; netcdf = True

#==============================================================================
# Fichiers c++
#==============================================================================
cpp_srcs =  ['Converter/Converter1.cpp',
             'Converter/copy.cpp',
             'Converter/extractVars.cpp',
             'Converter/initVars.cpp',
             'Converter/addVars.cpp',
             'Converter/randomizeVar.cpp',
             'Converter/Converter7.cpp',
             'Converter/ediff.cpp',
             'Converter/norm.cpp',
             'Converter/normalize.cpp',
             'Converter/magnitude.cpp',
             'Converter/isFinite.cpp',
             'Converter/setNANValuesAt.cpp',
             'Converter/convertBAR2Struct.cpp',
             'Converter/convertStruct2Tetra.cpp',
             'Converter/convertStruct2TetraBary.cpp',
             'Converter/convertStruct2Hexa.cpp',
             'Converter/convertStruct2NGon.cpp',
             'Converter/convertHexa2Struct.cpp',
             'Converter/convertUnstruct2NGon.cpp',
             'Converter/convertUnstruct2Hexa.cpp',
             'Converter/convertHexa2Tetra.cpp',
             'Converter/convertPrism2Tetra.cpp',
             'Converter/convertPyra2Tetra.cpp',
             'Converter/convertNGon2TetraBary.cpp',
             'Converter/convertArray2TetraBary.cpp',
             'Converter/convertHO2LO.cpp',
             'Converter/convertLO2HO.cpp',
             'Converter/convertTri2Quad.cpp',
             'Converter/convertQuad2Tri.cpp',
             'Converter/convertMix2BE.cpp',
             'Converter/convertStrand2Penta.cpp',
             'Converter/convertPenta2Strand.cpp',
             'Converter/center2Node.cpp',
             'Converter/center2Node_OLD.cpp',
             'Converter/node2Center.cpp',
             'Converter/node2Center_OLD.cpp',
             'Converter/node2ExtCenter.cpp',
             'Converter/extCenter2Node.cpp',
             'Converter/center2ExtCenter.cpp',
             'Converter/convertFilePyTree.cpp',
             'Converter/convertFilePyTreeTau.cpp',
             'Converter/convertFilePyTreeFsdm.cpp',
             'Converter/setPartialFields.cpp',
             'Converter/setPartialFieldsToSum.cpp',
             'Converter/filterPartialFields.cpp',
             'Converter/sendRecv.cpp',
             'Converter/IO/DynArrayIO.cpp',
             'Converter/IO/GenIO.cpp',
             'Converter/IO/GenIO_endian.cpp',
             'Converter/IO/GenIO_fmt.cpp',
             'Converter/IO/GenIO_fmttp.cpp',
             'Converter/IO/GenIO_fmtv3d.cpp',
             'Converter/IO/GenIO_fmtpov.cpp',
             'Converter/IO/GenIO_fmtmesh.cpp',
             'Converter/IO/GenIO_fmtgmsh.cpp',
             'Converter/IO/GenIO_bingmsh.cpp',
             'Converter/IO/GenIO_fmtobj.cpp',
             'Converter/IO/GenIO_binstl.cpp',
             'Converter/IO/GenIO_binply.cpp',
             'Converter/IO/GenIO_fmtstl.cpp',
             'Converter/IO/GenIO_fmtselig.cpp',
             'Converter/IO/GenIO_bin3ds.cpp',
             'Converter/IO/GenIO_bintp.cpp',
             'Converter/IO/GenIO_bince.cpp',
             'Converter/IO/GenIO_fmtplot3d.cpp',
             'Converter/IO/GenIO_binplot3d.cpp',
             'Converter/IO/GenIO_bintp108.cpp',
             'Converter/IO/GenIO_bintp75.cpp',
             'Converter/IO/GenIO_binv3d.cpp',
             'Converter/IO/GenIO_bindf3.cpp',
             'Converter/IO/GenIO_binwav.cpp',
             'Converter/IO/GenIO_binvtk.cpp',
             'Converter/IO/GenIO_fmtxfig.cpp',
             'Converter/IO/GenIO_fmtsvg.cpp',
             'Converter/IO/GenIO_fmtgts.cpp',
             'Converter/IO/GenIO_fmtcedre.cpp',
             'Converter/IO/GenIO_binarc.cpp',
             'Converter/IO/GenIO_fmtfoam.cpp',
             'Converter/IO/GenIO_fmtSU2.cpp',
             'Converter/IO/GenIO_bingltf.cpp',
             'Converter/IO/GenIO_createElts.cpp',
             'Converter/IO/GenIO_cplot.cpp',
             'Converter/IO/getBCFaces.cpp',
             'Converter/IO/convertPyTree2FFD.cpp',
             'Converter/cpyGhost2Real.cpp',
             'Converter/convertArray2Node.cpp',
             'Converter/detectEmptyBC.cpp',
             'Converter/tagDefinedBC.cpp',
             'Converter/fillJoin.cpp',
             'Converter/fillJoinNM.cpp',
             'Converter/fillCornerGhostCells.cpp',
             'Converter/fillCornerGhostCells2.cpp',
             'Converter/getBorderIndices.cpp',
             'Converter/getJoinDonorIndices.cpp',
             'Converter/conformizeNGon.cpp',
             'Converter/conformizeNGon1.cpp',
             'Converter/identify.cpp',
             'Converter/nearest.cpp',
             'Converter/identifySolutions.cpp',
             'Converter/hook.cpp',
             'Converter/globalHook.cpp',
             'Converter/globalIndex.cpp',
             'Converter/createBBTree.cpp',
             'Converter/ADF/ADF_interface.cpp',
             'Converter/ADF/ADF_internals.cpp',
             'Converter/ADF/cgns_io.cpp',
             'Converter/addGhostCellsNGon.cpp',
             'Converter/rmGhostCellsNGon.cpp',
             'Converter/extractBCMatch.cpp',
             'Converter/Adapter/adaptPE2NFace.cpp',
             'Converter/Adapter/adaptNFace2PE.cpp',
             'Converter/Adapter/adaptBCFace2BCC.cpp',
             'Converter/Adapter/adaptBCC2BCFace.cpp',
             'Converter/Adapter/adaptBCFacePL2VertexPL.cpp',
             'Converter/Adapter/adaptNGon2Index.cpp',
             'Converter/Adapter/adaptNFace2Index.cpp',
             'Converter/Adapter/adaptNGon42NGon3.cpp',
             'Converter/Adapter/adaptNGon32NGon4.cpp',
             'Converter/Adapter/adaptShiftedPE2PE.cpp',
             'Converter/Adapter/signNGonFaces.cpp',
             'Converter/Adapter/unsignNGonFaces.cpp',
             'Converter/Adapter/makeParentElements.cpp',
             'Converter/Adapter/adaptSurfaceNGon.cpp',
             'Converter/Adapter/adapt2FastP.cpp',
             'Converter/Adapter/createElsaHybrid.cpp',
             'Converter/Adapter/pointList2Ranges.cpp',
             'Converter/Adapter/pointList2SPL.cpp',
             'Converter/Adapter/range2PointList.cpp',
             'Converter/Adapter/PR2VL.cpp',
             'Converter/Adapter/diffIndex.cpp',
             'Converter/setBCDataInGhostCells.cpp',
             'Converter/extrapInterior2BCFace.cpp',
             'Converter/nullifyVectorAtBCFace.cpp',
             'Converter/nuga_ghost.cpp',
             'Converter/extractBCFields.cpp',
             'Converter/Extract/extractFields.cpp']
cpp_srcs += ['Converter/IO/GenIO_adfcgns.cpp']

if EXPRESSION:
    cpp_srcs += ['Converter/Expression/ast.cpp',
                 'Converter/Expression/function.cpp',
                 'Converter/Expression/lexer.cpp',
                 'Converter/Expression/math_function.cpp',
                 'Converter/Expression/parser.cpp',
                 'Converter/Expression/symbol_table.cpp',
                 'Converter/Expression/simd_vector_wrapper.cpp']

if hdf:
    cpp_srcs += ['Converter/IO/GenIO_hdfcgns.cpp',
                 'Converter/IO/GenIO_hdffsdm.cpp']
else:
    cpp_srcs += ['Converter/IO/GenIO_hdfcgns_stub.cpp',
                 'Converter/IO/GenIO_hdffsdm_stub.cpp']

if netcdf:
    cpp_srcs += ['Converter/IO/GenIO_tau.cpp']
else:
    cpp_srcs += ['Converter/IO/GenIO_tau_stub.cpp']

# png
cpp_srcs += ['Converter/IO/GenIO_binpng.cpp']
# jpg
cpp_srcs += ['Converter/IO/GenIO_binjpg.cpp']

#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ['Converter/Fortran/WriteBCFileF.for',
            'Converter/Fortran/WriteGridF.for',
            'Converter/Fortran/WriteGridHeadF.for',
            'Converter/Fortran/WriteGridF.for',
            'Converter/Fortran/WriteIBFileF.for',
            'Converter/Fortran/WriteFFDFileF.for']

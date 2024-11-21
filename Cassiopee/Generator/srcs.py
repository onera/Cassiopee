NETGEN = True; TETGEN = True; MMGS = True

#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ["Generator/cart.cpp",
            "Generator/cylinder.cpp",
            "Generator/generator3.cpp",
            "Generator/front2Hexa.cpp",
            "Generator/front2Struct.cpp",
            "Generator/hyper2d.cpp",
            "Generator/getNormalMap.cpp",
            "Generator/getVolumeMap.cpp",
            "Generator/getCellCenters.cpp",
            "Generator/getFaceCentersAndAreas.cpp",
            "Generator/getOrthogonalityMap.cpp",
            "Generator/getRegularityMap.cpp",
            "Generator/getAngleRegularityMap.cpp",
            "Generator/getCircumCircleMap.cpp",
            "Generator/getInCircleMap.cpp",
            "Generator/barycenter.cpp",
            "Generator/checkPointInCEBB.cpp",
            "Generator/enforce.cpp",
            #"Generator/tanhdist.cpp",
            "Generator/enforceCurvature.cpp",
            "Generator/checkMesh.cpp",
            "Generator/getCEBBIntersection.cpp",
            "Generator/bboxIntersection.cpp",
            "Generator/obboxIntersection.cpp",
            #"Generator/obboxTreeIntersection.cpp",
            "Generator/getBBoxOfCells.cpp",
            "Generator/obbox.cpp",
            "Generator/computeCellPlanarity.cpp",
            "Generator/map.cpp",
            "Generator/elliptic.cpp",
            "Generator/grow.cpp",
            "Generator/stack.cpp",
            "Generator/cartHexa.cpp",
            "Generator/cartTetra.cpp",
            "Generator/cartPenta.cpp",
            "Generator/cartNGon.cpp",
            "Generator/cartPyra.cpp",
            "Generator/cartr1.cpp",
            "Generator/cartr2.cpp",
            "Generator/close.cpp",
            "Generator/closeBorders.cpp",
            "Generator/pointedHat.cpp",
            "Generator/stitchedHat.cpp",
            "Generator/delaunay.cpp",
            "Generator/checkDelaunay.cpp",
            "Generator/TFI.cpp",
            "Generator/TFI2D.cpp",
            "Generator/TFI3D.cpp",
            "Generator/TFITRI.cpp",
            "Generator/TFIPENTA.cpp",
            "Generator/TFITETRA.cpp",
            "Generator/selectInsideElts.cpp",
            "Generator/densifyMesh.cpp",
            "Generator/getTriQualityMap.cpp",
            "Generator/T3mesher2D.cpp",
            "Generator/fittingPlaster.cpp",
            "Generator/gapfixer.cpp",
            "Generator/gapsmanager.cpp",
            "Generator/octree.cpp",
            "Generator/octree3.cpp",
            "Generator/octree2AMR.cpp",
            "Generator/extendCartGrids.cpp",
            "Generator/octree2Struct.cpp",
            "Generator/balanceOctree.cpp",
            "Generator/adaptOctree.cpp",
            "Generator/adaptOctree3.cpp",
            "Generator/conformOctree3.cpp",
            "Generator/expandLayer.cpp",
            "Generator/octreeTbx.cpp",
            "Generator/snapFront.cpp",
            "Generator/snapSharpEdges.cpp",
            "Generator/surfaceWalk.cpp",
            "Generator/extrusionTbx.cpp",
            "Generator/getEdgeRatio.cpp",
            "Generator/getMaxLength.cpp",
            "Generator/quad2pyra.cpp",
            "Generator/blankSelf.cpp",
            "Generator/extrapWithCellN.cpp"
            ]

cpp_srcs2 = []
# netgen
if NETGEN:
    srcs1 = ["adfront2.cpp", "adfront3.cpp", "bisect.cpp", "boundarylayer.cpp", 
             "clusters.cpp", "curvedelems.cpp", "delaunay.cpp", 
             "delaunay2d.cpp", "geomsearch.cpp", "global.cpp", 
             "hprefinement.cpp", "improve2.cpp",
             "improve2gen.cpp", "improve3.cpp", "localh.cpp", "meshclass.cpp",  
             "meshfunc.cpp", "meshfunc2d.cpp", "meshing2.cpp", "meshing3.cpp",  
             "meshtool.cpp", "meshtype.cpp", "msghandler.cpp", "netrule2.cpp",  
             "netrule3.cpp", "parser2.cpp", "parser3.cpp", "prism2rls.cpp",  
             "pyramid2rls.cpp", "pyramidrls.cpp", "quadrls.cpp", "refine.cpp", 
             "ruler2.cpp", "ruler3.cpp", "secondorder.cpp", 
             "smoothing2.5.cpp",  
             "smoothing2.cpp", "smoothing3.cpp", "specials.cpp", 
             "tetrarls.cpp",  
             "topology.cpp", "triarls.cpp", "validate.cpp", "zrefine.cpp", 
             "bcfunctions.cpp",
             "parallelmesh.cpp", "paralleltop.cpp",
             "basegeom.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/meshing/"+i]
    srcs1 = ["array.cpp", "bitarray.cpp", "dynamicmem.cpp", "flags.cpp",
             "hashtabl.cpp", "mystring.cpp", "ngexception.cpp", "optmem.cpp",
             "parthreads.cpp", "profiler.cpp", "seti.cpp", "sort.cpp", 
             "spbita2d.cpp", "symbolta.cpp", "table.cpp",
             "mpi_interface.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/general/"+i]
    srcs1 = ["densemat.cpp", "polynomial.cpp",  "bfgs.cpp", 
             "linopt.cpp", "linsearch.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/linalg/"+i]
    srcs1 = ["adtree.cpp", "geom2d.cpp", "geom3d.cpp", "geomfuncs.cpp",
             "geomtest3d.cpp", "transform3d.cpp", "spline.cpp", 
             "splinegeometry.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/gprim/"+i]
    srcs1 = ["genmesh2d.cpp", "geom2dmesh.cpp", "geometry2d.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/geom2d/"+i]
    srcs1 = ["meshstlsurface.cpp", "stlgeom.cpp", "stlline.cpp",
             "stltool.cpp", "stlgeomchart.cpp", "stlgeommesh.cpp",
             "stltopology.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/stlgeom/"+i]
    srcs1 = ["algprim.cpp", "brick.cpp", "bspline2d.cpp", "csgeom.cpp",
             "csgparser.cpp", "curve2d.cpp", "edgeflw.cpp",
             "explicitcurve2d.cpp", "extrusion.cpp", "gencyl.cpp",
             "genmesh.cpp", "identify.cpp", "manifold.cpp", "meshsurf.cpp",
             "polyhedra.cpp", "revolution.cpp", "singularref.cpp",
             "solid.cpp", "specpoin.cpp", "spline3d.cpp", "surface.cpp",
             "triapprox.cpp"]
    for i in srcs1: cpp_srcs2 += ["Generator/Netgen/csg/"+i]

    cpp_srcs2 += ["Generator/netgen.cpp"]
else:
    cpp_srcs2 += ["Generator/netgen_stub.cpp"]

# tegen
if TETGEN:
    cpp_srcs2 += ["Generator/Tetgen/tetgen.cxx",
                  "Generator/Tetgen/predicates.cxx",
                  "Generator/tetgen.cpp"]
else:
    cpp_srcs2 += ["Generator/tetgen_stub.cpp"]

# mmg
if MMGS:
    cpp_srcs2 += ["Generator/mmgs.cpp",
                  "Generator/MMGS/mmg.c",    
                  "Generator/MMGS/mmgs1.c",
                  "Generator/MMGS/mmgs2.c",
                  "Generator/MMGS/analys_s.c",
                  "Generator/MMGS/anisomovpt.c",
                  "Generator/MMGS/anisomovpt_s.c",
                  "Generator/MMGS/anisosiz.c",
                  "Generator/MMGS/anisosiz_s.c",
                  "Generator/MMGS/API_functions.c",
                  "Generator/MMGS/API_functions_s.c",
                  "Generator/MMGS/bezier.c",
                  "Generator/MMGS/bezier_s.c",
                  "Generator/MMGS/boulep.c",
                  "Generator/MMGS/boulep_s.c",
                  "Generator/MMGS/chkmsh_s.c",
                  "Generator/MMGS/chrono.c",
                  "Generator/MMGS/colver_s.c",
                  "Generator/MMGS/eigenv.c",
                  "Generator/MMGS/gentools_s.c",
                  "Generator/MMGS/hash.c",
                  "Generator/MMGS/hash_s.c",
                  "Generator/MMGS/inout.c",
                  "Generator/MMGS/inout_s.c",
                  "Generator/MMGS/intmet.c",
                  "Generator/MMGS/intmet_s.c",
                  "Generator/MMGS/isosiz.c",
                  "Generator/MMGS/isosiz_s.c",
                  "Generator/MMGS/libmmgs.c",
                  "Generator/MMGS/libmmgs_tools.c",
                  "Generator/MMGS/librnbg.c",
                  "Generator/MMGS/librnbg_s.c",
                  "Generator/MMGS/mettools.c",
                  "Generator/MMGS/movpt_s.c",
                  "Generator/MMGS/quality.c",
                  "Generator/MMGS/quality_s.c",
                  "Generator/MMGS/scalem.c",
                  "Generator/MMGS/split_s.c",
                  "Generator/MMGS/swapar_s.c",
                  "Generator/MMGS/tools.c",
                  "Generator/MMGS/variadic_s.c",
                  "Generator/MMGS/zaldy_s.c",
                  "Generator/MMGS/mmgexterns.c",
                  "Generator/MMGS/mmgsexterns.c"]
else:
    cpp_srcs2 += ["Generator/mmgs_stub.cpp"]

#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ['Generator/Fortran/MeshqualF.for',
            'Generator/Fortran/HgF.for',
            'Generator/Fortran/Hg2F.for',
            'Generator/Fortran/Hg3F.for',
            'Generator/Fortran/Hg4F.for',
            'Generator/Fortran/DecbtF.for',
            'Generator/Fortran/SolbtF.for',
            'Generator/Fortran/PtridF.for',
            'Generator/Fortran/DecF.for',
            'Generator/Fortran/SolF.for',
            'Generator/Fortran/GelgF.for',
            'Generator/Fortran/StretchF.for',
            'Generator/Fortran/MetsF.for',
            'Generator/Fortran/CoeffttmF.for',
            'Generator/Fortran/Sor9ptF.for']

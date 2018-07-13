#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ["Geom/naca.cpp",
            "Geom/Geom1.cpp",
            "Geom/Geom2.cpp",
            "Geom/lineGenerate.cpp",
            "Geom/lineGenerate2.cpp",
            "Geom/axisym.cpp",
            "Geom/getCurvatureAngle.cpp",
            "Geom/getCurvatureRadius.cpp",
            "Geom/getCurvatureHeight.cpp",
            "Geom/polyline.cpp",
            "Geom/spline.cpp",
            "Geom/nurbs.cpp",
            "Geom/bezier.cpp",
            "Geom/volumeFromCrossSections.cpp",
            "Geom/addSeparationLine.cpp",
            "Geom/torus.cpp",
            "Geom/getSharpestAngleForVertices.cpp",
            "Geom/getNearestPointIndex.cpp"]

#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ['Geom/Fortran/NacaF.for',
            'Geom/Fortran/nacaS4GeneF.for',
            'Geom/Fortran/nacaS5GeneF.for',
            'Geom/Fortran/nacaS4ModGeneF.for',
            'Geom/Fortran/AxisymF.for']

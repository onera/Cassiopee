// This is a povray sample file
//
// Please render this file with :
// povray -W800 -H600 +a0.3 +SP16 render.pov +P
//=============================================================================
// General includes
//=============================================================================
#include "colors.inc"
#declare Blue_1   = rgb <0.263,0.345,0.549>;
#declare Blue_2   = rgb <0.243,0.451,0.765>;
#declare Red_1    = rgb <0.761,0.278,0.220>;
#declare Orange_1 = rgb <0.776,0.482,0.212>;
#declare Yellow_1 = rgb <0.757,0.604,0.231>;
#declare Green_1  = rgb <0.282,0.455,0.216>;
#declare Green_2   = rgb <0.675,0.671,0.416>;
#declare Rose_1   = rgb <0.694,0.549,0.514>;
#declare Rose_2   = rgb <0.871,0.702,0.675>;
//=============================================================================
// Include the mesh
//=============================================================================
#include "out.pov"
object {out 
        translate -0.3*z
          finish { phong 0.9  phong_size 60 metallic }  
          pigment {Rose_2}
}
//=============================================================================
// Camera position
//=============================================================================
camera {
  perspective
    location <-1.,1.,0.2>
    look_at <0.4,0.1,-0.1>
    sky <0,0,1>
    }
//=============================================================================
// Sol
//=============================================================================
//plane {z,-1 pigment{Green_2}}
plane {z,-0.5
         texture {
  pigment {color rgbf < 0.1, 0.1, 1, 1>}
  finish {
    ambient 0 diffuse 0.7
      reflection {0.3, 1
                    metallic //use metallic reflection
                    }
    conserve_energy
      metallic //use metallic highlights
      }
  normal {bumps bump_size 0.175 scale < 4, 1, 1>*0.025}
}
}
//=============================================================================
// Sky
//=============================================================================
//text {
//    ttf "timrom.ttf" "Cassiopee" 0.25, 0
//      finish { reflection .25 specular 1 }
//    rotate <0,180,-180.>
//      translate <0,0,-0.8>
//      
//    pigment { White }
//  }
//sky_sphere {pigment {rgb <0.0,0.910,0.831>}}
sky_sphere {pigment {rgb <0.0,0.610,0.831>}}
light_source {<-15,-5,5> White*1. parallel point_at <0,0,0.4>}
light_source {<5,5,2> 2*White}
//light_source {<-16,-17,16> Green*2}

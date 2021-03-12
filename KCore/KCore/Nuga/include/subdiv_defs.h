/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef SUBDIV_DEFS_H
#define SUBDIV_DEFS_H

namespace NUGA
{
  enum eSUBDIV_TYPE { ISO = 0, ISO_HEX=1/*all to hexa*/, ISO2/*iso metric field : spheres*/, DIR, ANISO/*aniso medtric field : ellipses*/ };
  enum eDIR { NONE = 0, Xd/*d because conflict with eClassifyer::X*/, Y, XY, /*XZ, YZ*/XYZ };

  typedef void(*reordering_func)(E_Int* child, bool reverse, E_Int i0);
}

#endif

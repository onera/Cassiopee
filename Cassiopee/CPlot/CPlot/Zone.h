/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _CPLOT_ZONE_H_
#define _CPLOT_ZONE_H_
#include "ZoneImpl.h"
#include <vector>

#ifndef MAXSTRINGLENGTH
#define MAXSTRINGLENGTH 128   /* max length of strings */
#endif

/* Define a generic zone.
 Can be derived in structured and unstructured zones. */
class Zone
{
  public:
    // fast inverse root square
    inline float invsqrt( float x ) 
    {
        union intfloat { int   i; float x; };
        intfloat var;
        var.x = x;
        float xhalf = 0.5f * x;
        var.i = 0x5f375a86 - ( var.i >> 1 ); // gives initial guess y0
        var.x = var.x * ( 1.5f - xhalf * var.x * var.x ); // Newton step, repeating increases accuracy
        return var.x;
    };

    // Steal the implementation pointer !
    Zone(CPlotState* states, ZoneImpl* impl);
    virtual ~Zone();

    CPlotState* ptrState;
    char zoneName[MAXSTRINGLENGTH];  // zone name
    
    E_Int dim;                       // dim = 1, 2 or 3
    E_Int nfield;                    // Number of variables defined in zone
                                     // that is 3 (xyz) + nfield
    E_Int npts;                      // number of vertices
    double* x;                       // storage of coords of points (x,y,z)
    double* y;
    double* z;
    double **f;                      // storage of variables (size nf)
    char** varnames;                 // x,y,z + other variables
    std::vector<float*> surf;        // normal to each element
                                     // for struct: only one at faces
                                     // for TRI, QUAD: only one at nodes
                                     // for others, one for each connectivity at faces
    double xmin, xmax;               // min max of zone coords
    double ymin, ymax; 
    double zmin, zmax;
    double* minf;                    // min max value of f (size nf)
    double* maxf;
    double Di, Dj, Dk;                // Mesh densities
    double xc, yc, zc;                // center of zone (sphere)
    double rc;                        // radius of sphere
    E_Int active; // if set to 0 the zone is not displayed
    E_Int selected; // if set to 1 the zone is displayed selected
    E_Int previouslySelected; // 1 means that zone was selected at previous click
    E_Int blank; // no blanking (0), otherwise the blanking variable number

    char renderTag[MAXSTRINGLENGTH];  // tag for render
    float colorR, colorG, colorB; // color RGB (from tag)
    E_Int material; // material code (from tag)
    double blending; // blending (from tag)
    E_Int meshOverlay; // 1 if mesh must be overlaid (from tag)
    float meshColorR, meshColorG, meshColorB; // mesh color (from tag)
    float meshWidth; // width of mesh lines (from tag)
    float shaderParam1, shaderParam2; // shaders param (from tag)
                                  // default value is 1. (can vary in 0-2).
    double* texu; // pointe sur le champ u pour les textures
    double* texv; // pointe sur le champ v pour les textures
    double* texw; // pointe sur le champ w pour les textures
    double* regtexu; // generated regtexu if not in fields
    double* regtexv; // generated regtexv if not in fields
    float tesOuter, tesInner; // Tesselation shader param ( from tag ? )
                              // Default value is 3 ( can vary in 3-50 )

    ZoneImpl* ptr_impl;
    void setUseGPURessources()   { ptr_impl->_GPUResUse = 1; }
    void unsetUseGPURessources() { ptr_impl->_GPUResUse = 0; }
    void freeGPURessources(bool useGPURessources=true, bool freeIso=true);
    void destroyGPUIsoField() { ptr_impl->destroyIsoField(); }
    unsigned char* _voxelArray; // array voxel pour les zones en rendu volumique
};

#endif

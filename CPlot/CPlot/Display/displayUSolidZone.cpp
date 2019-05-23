/*    
    Copyright 2013-2019 Onera.

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
#include "../Data.h"

//=============================================================================
/*
  Display une zone en solid ou en material (mode SOLID ou RENDER)
*/
//=============================================================================
void Data::displayUSolidZone( UnstructZone *zonep, int zone, int zonet )
{
    int i, n1, n2, n3, n4, n5, n6, n7, n8;
    int ret1, ret2, ret3, ret4, f;

    // Style
    float color1[3]; float color2[ 3 ];

    // Colormap
    float r, g, b;
    void (*getrgb)(Data* data, double, float*, float*, float*);
    getrgb = _plugins.colorMap->next->f;

    E_Float nz = 1. / _numberOfUnstructZones;
#include "solidStyles.h"

    // Ecrasement si renderTag
    if ( zonep->colorR > -0.5 ) {
        color1[ 0 ] = zonep->colorR;
        color1[ 1 ] = zonep->colorG;
        color1[ 2 ] = zonep->colorB;
    }

#include "selection.h"

    bool is1D = ( ( zonep->eltType == 1 ) | ( zonep->eltType == 10 && zonep->nelts1D > 0 ) );
    if (is1D == true && ptrState->mode == RENDER)
        glLineWidth( 1. + 5 * zonep->shaderParam1 );
    else if (is1D == true) glLineWidth(3.);
    else glLineWidth(1.);

    // scale
    E_Float s = MAX( zonep->xmax - zonep->xmin, zonep->ymax - zonep->ymin );
    s = MAX( s, zonep->zmax - zonep->zmin );
    s = 100. / ( s + 1.e-12 );

    // Only for textured rendering, we use vect display =======================
    if (ptrState->mode == RENDER && zonep->material == 14 && zonep->texu != NULL)  // Textured rendering
    {
#ifdef __SHADERS__
        triggerShader(*zonep, zonep->material, s, color1);
#endif

        int ff;
        double offb = 0.;
#undef PLOT
        double* f1 = zonep->texu;
        double* f2 = zonep->texv;
        double* f3 = zonep->texw;
        double fmin1, fmax1, fmin2, fmax2, fmin3, fmax3;
        fmax1 = 0.; fmin1 = 1.;
        fmax2 = 0.; fmin2 = 1.;
        fmax3 = 0.; fmin3 = 1.;
#define GL_QUADS_ARE GL_QUADS
#define PLOTQUAD PLOTQUADQ
#define PLOTQUAD2 PLOTQUADQ2
#include "displayUVectSolidZone.h"
        glLineWidth( 1. );
        return;
    }
    // END Textured rendering ============================================
#ifdef __SHADERS__
    if ( ptrState->mode == RENDER ) {
        if ( zonep->selected == 1 && zonep->active == 1 )
            triggerShader( *zonep, zonep->material, s, color2 );
        else
            triggerShader( *zonep, zonep->material, s, color1 );
    } else {
        if ( zonep->selected == 1 && zonep->active == 1 )
            triggerShader( *zonep, 0, s, color2 );
        else
            triggerShader( *zonep, 0, s, color1 );
    }
#endif

#undef PLOT
#undef PLOTTRI
#undef PLOTTRI2
#undef PLOTQUAD
#undef PLOTQUAD2
#include "displayUSolidZone.h"
    glLineWidth( 1. );
}

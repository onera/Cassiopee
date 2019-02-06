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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//=============================================================================
/* 
   Display unstructured meshes (mode SOLID ou RENDER)
   This surface mesh is plotted with plane quads or tri.
*/
//=============================================================================
void DataDL::displayUSolid()
{
    int zone, zonet;
    if ( _numberOfUnstructZones == 0 ) return;

    // Enable blending
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glEnable( GL_POLYGON_OFFSET_FILL );
    glEnable( GL_POLYGON_OFFSET_LINE );

    zone = 0;
    while ( zone < _numberOfUnstructZones ) {
        glPolygonOffset( 1., zone % 10 + 1. );
        zonet = zone + _numberOfStructZones;
        UnstructZone *zonep = _uzones[ zone ];
        ZoneImplDL *zoneImpl = static_cast<ZoneImplDL *>( zonep->ptr_impl );

        // if zone is active and in frustum
        if ( ( zonep->active == 1 ||
               ( zonep->active == 0 && ptrState->ghostifyDeactivatedZones == 1 ) ) &&
             isInFrustum( zonep, _view ) == 1 ) {
            if ( ptrState->mode == RENDER && zonep->meshOverlay == 1 ) {
                noLight();
                _shaders.activate( (short unsigned int)0 );
                displayUMeshZone( zonep, zone, zonet );
                light( 2 );
            }

            double alphaSav = ptrState->alpha;
            if ( ptrState->mode == RENDER && zonep->blending != -1. ) {
                ptrState->alpha = zonep->blending;
                if ( ptrState->alpha < 0.9999 ) { /*glEnable(GL_CULL_FACE);*/
                    glDepthMask( GL_FALSE );
                }
            }

            if ( ptrState->mode == RENDER && zonep->colorR < -1.5 )  // Iso
            {
#ifdef __SHADERS__
#include "isoShaders.h"
#endif
                if ( ptrState->isoLight == 1 && ptrState->dim == 3 ) light( 3 );

                //displayUIsoSolidZone(zonep, zone, (int)(-zonep->colorR-2));
                renderUIsoSolidZone( zonep, zone, (int)( -zonep->colorR - 2 ) );
            }
            //else if (ptrState->mode == RENDER && zonep->material == 9) // Billboarding
            //{
            //displayBillBoards((Zone*)zonep, zone);
            //}
            else  // Autres cas
            {
                if ( ptrState->dim == 3 ) {
                    switch ( ptrState->solidStyle ) {
                        case 0:  // one side
                            light( 0 );
                            break;

                        case 1:  // two sides
                            light( 2 );
                            break;

                        case 3:  // two sides
                            light( 2 );
                            break;

                        case 4:  // two sides
                            light( 2 );
                            break;
                    }
                }
                if ( zonep->_is_high_order == false ) {
                    if ( zoneImpl->_DLsolid != 0 )
                        renderGPUUSolidZone( zonep, zone, zonet );
                    else
                        displayUSolidZone( zonep, zone, zonet );  // Direct
                } else {
                    if ( zoneImpl->_DLsolid != 0 )
                        renderGPUUSolidHOZone( zonep, zone, zonet );
                    else
                        displayUSolidHOZone( zonep, zone, zonet );  // Direct
                }
            }
            ptrState->alpha = alphaSav;
            glDisable( GL_CULL_FACE );
            glDepthMask( GL_TRUE );
        }
        zone++;
    }

    noLight();
#ifdef __SHADERS__
    _shaders.activate( (short unsigned int)0 );
    glActiveTexture( GL_TEXTURE1 );
    glDisable( GL_TEXTURE_2D );
    glActiveTexture( GL_TEXTURE0 );
    glDisable( GL_TEXTURE_2D );
#endif
    glDisable( GL_BLEND );
    glDisable( GL_POLYGON_OFFSET_FILL );
    glDisable( GL_POLYGON_OFFSET_LINE );
    glColor4f( 1., 1., 1., 1. );
}

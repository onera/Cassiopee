/*    
    Copyright 2013-2026 ONERA.

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

#ifndef _CPLOT_PREFERENCES_H_
#define _CPLOT_PREFERENCES_H_

/* Preferences */
typedef struct
{
    /* -- Performance preferences -- */
    // Speed factor (0 is accurate drawing)
    E_Int speed;
    // Number of roll steps when rotating views
    E_Int nroll;

    /* -- Smoke preferences -- */
    // Max number of particles in smoke
    E_Int maxParticles;
    // Radius of smoke particles
    float smokeRadius;
    // Radius of emission zone
    float emissionRadius; 
    // Time step for smoke
    float smokeTimeStep; 

    /* -- Texture preferences -- */
    // Number of mesh quads to map texture
    E_Int sizeTexi;
    E_Int sizeTexj;

    /* -- Plugins preferences -- */
    // Plugins location
    char pluginsLocation[MAXSTRINGLENGTH];
    // Current colormap function
    struct chain_function_double* colorMap;
    // Current look for function
    struct chain_function_void* lookFor;
    // Current add-a-variable function
    struct chain_function_void* addAVariable;
    // Current blanking function
    struct chain_function_int* blanking;
    // Current selection function
    struct chain_function_void* selectNextZone;
    struct chain_function_void* selectPreviousZone;
    // Current mouse function
    struct chain_function_void3* mouseClick;
    struct chain_function_void3* mouseMultipleClick;
    struct chain_function_void3* mouseAccurateClick;
    struct chain_function_void3* mouseRightClick;
    // Current screen dump function
    struct chain_function_void2* screenDump;

} Preferences;

// Load prefs
void loadPrefs();

#endif

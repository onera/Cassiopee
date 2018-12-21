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

#ifndef _CPLOT_PARTICLES_H_
#define _CPLOT_PARTICLES_H_

typedef struct
{
    /* Center of particle */
    float x, y, z;
    /* Concentration of smoke in particle */
    float conc;
    float cinit;
    /* Color of particule */
    float r, g, b;
    /* Previous displacement */
    float dx, dy, dz;
    /* Radius */
    float rad;
    /* Current interpolation domain */
    int interpDomain; 
} Particle;

typedef struct 
{
    /* Initial position of smoke emission */
    float x, y, z;
    /* Start domain */
    int startDomain;
    /* Number of particles */
    int n;
    /* Particles */
    Particle* particle;
    /* Variables location */
    int pro;
    int prou;
    int prov;
    int prow;
    int pcellN;
    /* Emission zone radius */
    float de;
} Cloud;

// Called by cplot for creating smoke
void initCloud();
void launchParticles();
void initParticles();
void initParticle(Particle *p);
void drawParticles(float minlife);
int setCheckTexture();

#endif

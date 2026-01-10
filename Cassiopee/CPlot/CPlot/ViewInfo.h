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

#ifndef _CPLOT_VIEWINFO_H_
#define _CPLOT_VIEWINFO_H_

#include "cplot.h"

#define ANG2RAD 0.017453292519943295
#define PIS2 1.5707963267948966

/* Define informations on view */
typedef struct
{
    double xcam;   // camera position
    double ycam;
    double zcam;
    double xeye;   // position to where is pointing the eye 
    double yeye;
    double zeye;
    double dirx;   // camera direction
    double diry;
    double dirz;
    double xcam3D; // save the 3D camera position if set
    double ycam3D; 
    double zcam3D;
    double xeye3D; // save the 3D eye position if set
    double yeye3D; 
    double zeye3D;
    double dirx3D; // save the 3D camera direction if set
    double diry3D;
    double dirz3D;
    double xcam2D; // save the 2D position if set
    double ycam2D;
    double zcam2D;
    double xeye2D; // save the 2D eye position if set
    double yeye2D; 
    double zeye2D;
    double dirx2D; // save the 2D camera direction if set
    double diry2D;
    double dirz2D;
    double xcam1D; // save the 1D position if set
    double ycam1D;
    double zcam1D;
    double xeye1D; // save the 1D eye position if set
    double yeye1D; 
    double zeye1D;
    double dirx1D; // save the 1D camera direction if set
    double diry1D;
    double dirz1D;
    E_Int w;         // window size
    E_Int h;
    E_Int wSav;      // Store window size when in full screen
    E_Int hSav;

    // frustum info
    E_Int clipping; // 0: far, 1: near, 2: very close
    double ratio;   
    double angle;
    double nearD;
    double farD;
    double tang;
    double nh;
    double nw;
    double fh;
    double fw;
    
    // Near plane
    double NearD;
    double NearNx; double NearNy; double NearNz;
    // Far plane
    double FarD;
    double FarNx; double FarNy; double FarNz;
    // Top plane
    double TopD;
    double TopNx; double TopNy; double TopNz;
    // Bottom plane
    double BottomD;
    double BottomNx; double BottomNy; double BottomNz;
    // Left plane
    double LeftD;
    double LeftNx; double LeftNy; double LeftNz;
    // Right plane
    double RightD;
    double RightNx; double RightNy; double RightNz;
    
} ViewInfo;

#endif

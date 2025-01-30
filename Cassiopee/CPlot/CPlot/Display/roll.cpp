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
#include "../Data.h"

//=============================================================================
/* 
   Automatic movement to reach 2D saved camera position
   Current position is saved to 3D saved camera position
*/
//=============================================================================
void Data::roll3Dto2D()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz;
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Save the 3D camera position
  _view.xcam3D = _view.xcam;
  _view.ycam3D = _view.ycam;
  _view.zcam3D = _view.zcam;
  _view.xeye3D = _view.xeye;
  _view.yeye3D = _view.yeye;
  _view.zeye3D = _view.zeye;
  _view.dirx3D = _view.dirx;
  _view.diry3D = _view.diry;
  _view.dirz3D = _view.dirz;

  // Compute distance
  distx = _view.xcam - _view.xcam2D;
  disty = _view.ycam - _view.ycam2D;
  distz = _view.zcam - _view.zcam2D;
  edistx = _view.xeye - _view.xeye2D;
  edisty = _view.yeye - _view.yeye2D;
  edistz = _view.zeye - _view.zeye2D;
  ddistx = _view.dirx - _view.dirx2D;
  ddisty = _view.diry - _view.diry2D;
  ddistz = _view.dirz - _view.dirz2D;

  distx = distx*froll;
  disty = disty*froll;
  distz = distz*froll;
  edistx = edistx*froll;
  edisty = edisty*froll;
  edistz = edistz*froll;
  ddistx = ddistx*froll;
  ddisty = ddisty*froll;
  ddistz = ddistz*froll;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx;
    _view.ycam = _view.ycam - disty;
    _view.zcam = _view.zcam - distz;

    _view.xeye = _view.xeye - edistx;
    _view.yeye = _view.yeye - edisty;
    _view.zeye = _view.zeye - edistz;

    _view.dirx = _view.dirx - ddistx;
    _view.diry = _view.diry - ddisty;
    _view.dirz = _view.dirz - ddistz;    

    display();
  }
}

//=============================================================================
/*
  Automatic movement to reach 1D saved camera position
  Current position is saved to 2D saved camera position
*/
//=============================================================================
void Data::roll2Dto1D()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Save the 3D camera position
  _view.xcam2D = _view.xcam;
  _view.ycam2D = _view.ycam;
  _view.zcam2D = _view.zcam;
  _view.xeye2D = _view.xeye;
  _view.yeye2D = _view.yeye;
  _view.zeye2D = _view.zeye;
  _view.dirx2D = _view.dirx;
  _view.diry2D = _view.diry;
  _view.dirz2D = _view.dirz;

  // Compute distance
  distx = _view.xcam-_view.xcam1D;
  disty = _view.ycam-_view.ycam1D;
  distz = _view.zcam-_view.zcam1D;  
  edistx = _view.xeye-_view.xeye1D;
  edisty = _view.yeye-_view.yeye1D;
  edistz = _view.zeye-_view.zeye1D;
  ddistx = _view.dirx-_view.dirx1D;
  ddisty = _view.diry-_view.diry1D;
  ddistz = _view.dirz-_view.dirz1D;

  distx = distx*froll;
  disty = disty*froll;
  distz = distz*froll;
  edistx = edistx*froll;
  edisty = edisty*froll;
  edistz = edistz*froll;
  ddistx = ddistx*froll;
  ddisty = ddisty*froll;
  ddistz = ddistz*froll;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx;
    _view.ycam = _view.ycam - disty;
    _view.zcam = _view.zcam - distz;

    _view.xeye = _view.xeye - edistx;
    _view.yeye = _view.yeye - edisty;
    _view.zeye = _view.zeye - edistz;

    _view.dirx = _view.dirx - ddistx;
    _view.diry = _view.diry - ddisty;
    _view.dirz = _view.dirz - ddistz;    

    display();
  }
}

//=============================================================================
/*
  Automatic movement to reach 3D saved camera position
  Current position is saved to 1D saved camera position
*/
//=============================================================================
void Data::roll1Dto3D()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Save the 1D camera position
  _view.xcam1D = _view.xcam;
  _view.ycam1D = _view.ycam;
  _view.zcam1D = _view.zcam;
  _view.xeye1D = _view.xeye;
  _view.yeye1D = _view.yeye;
  _view.zeye1D = _view.zeye;
  _view.dirx1D = _view.dirx;
  _view.diry1D = _view.diry;
  _view.dirz1D = _view.dirz;

  // Compute distance
  distx = _view.xcam-_view.xcam3D;
  disty = _view.ycam-_view.ycam3D;
  distz = _view.zcam-_view.zcam3D;  
  edistx = _view.xeye-_view.xeye3D;
  edisty = _view.yeye-_view.yeye3D;
  edistz = _view.zeye-_view.zeye3D;
  ddistx = _view.dirx-_view.dirx3D;
  ddisty = _view.diry-_view.diry3D;
  ddistz = _view.dirz-_view.dirz3D;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx*froll;
    _view.ycam = _view.ycam - disty*froll;
    _view.zcam = _view.zcam - distz*froll;

    _view.xeye = _view.xeye - edistx*froll;
    _view.yeye = _view.yeye - edisty*froll;
    _view.zeye = _view.zeye - edistz*froll;

    _view.dirx = _view.dirx - ddistx*froll;
    _view.diry = _view.diry - ddisty*froll;
    _view.dirz = _view.dirz - ddistz*froll;    

    display();
  }
}

//=============================================================================
/*
  Automatic movement to reach 2D saved camera position
  Current position is saved to 1D saved camera position
*/
//=============================================================================
void Data::roll1Dto2D()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Save the 1D camera position
  _view.xcam1D = _view.xcam;
  _view.ycam1D = _view.ycam;
  _view.zcam1D = _view.zcam;
  _view.xeye1D = _view.xeye;
  _view.yeye1D = _view.yeye;
  _view.zeye1D = _view.zeye;
  _view.dirx1D = _view.dirx;
  _view.diry1D = _view.diry;
  _view.dirz1D = _view.dirz;

  // Compute distance
  distx = _view.xcam-_view.xcam2D;
  disty = _view.ycam-_view.ycam2D;
  distz = _view.zcam-_view.zcam2D;  
  edistx = _view.xeye-_view.xeye2D;
  edisty = _view.yeye-_view.yeye2D;
  edistz = _view.zeye-_view.zeye2D;
  ddistx = _view.dirx-_view.dirx2D;
  ddisty = _view.diry-_view.diry2D;
  ddistz = _view.dirz-_view.dirz2D;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx*froll;
    _view.ycam = _view.ycam - disty*froll;
    _view.zcam = _view.zcam - distz*froll;

    _view.xeye = _view.xeye - edistx*froll;
    _view.yeye = _view.yeye - edisty*froll;
    _view.zeye = _view.zeye - edistz*froll;

    _view.dirx = _view.dirx - ddistx*froll;
    _view.diry = _view.diry - ddisty*froll;
    _view.dirz = _view.dirz - ddistz*froll;    

    display();
  }
}

//=============================================================================
/*
  Automatic movement to reach 3D saved camera position
  Current position is saved to 2D saved camera position
*/
//=============================================================================
void Data::roll2Dto3D()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Save the 2D camera position
  _view.xcam2D = _view.xcam;
  _view.ycam2D = _view.ycam;
  _view.zcam2D = _view.zcam;
  _view.xeye2D = _view.xeye;
  _view.yeye2D = _view.yeye;
  _view.zeye2D = _view.zeye;
  _view.dirx2D = _view.dirx;
  _view.diry2D = _view.diry;
  _view.dirz2D = _view.dirz;

  // Compute distance
  distx = _view.xcam-_view.xcam3D;
  disty = _view.ycam-_view.ycam3D;
  distz = _view.zcam-_view.zcam3D;
  edistx = _view.xeye-_view.xeye3D;
  edisty = _view.yeye-_view.yeye3D;
  edistz = _view.zeye-_view.zeye3D;
  ddistx = _view.dirx-_view.dirx3D;
  ddisty = _view.diry-_view.diry3D;
  ddistz = _view.dirz-_view.dirz3D;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx*froll;
    _view.ycam = _view.ycam - disty*froll;
    _view.zcam = _view.zcam - distz*froll;
    
    _view.xeye = _view.xeye - edistx*froll;
    _view.yeye = _view.yeye - edisty*froll;
    _view.zeye = _view.zeye - edistz*froll;

    _view.dirx = _view.dirx - ddistx*froll;
    _view.diry = _view.diry - ddisty*froll;
    _view.dirz = _view.dirz - ddistz*froll;

    display();
  }
}

//=============================================================================
/* 
   Automatic movement to reach 1D saved camera position
   Current position is saved to 3D saved camera position
*/
//=============================================================================
void Data::roll3Dto1D()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Save the 3D camera position
  _view.xcam3D = _view.xcam;
  _view.ycam3D = _view.ycam;
  _view.zcam3D = _view.zcam;
  _view.xeye3D = _view.xeye;
  _view.yeye3D = _view.yeye;
  _view.zeye3D = _view.zeye;
  _view.dirx3D = _view.dirx;
  _view.diry3D = _view.diry;
  _view.dirz3D = _view.dirz;

  // Compute distance
  distx = _view.xcam-_view.xcam1D;
  disty = _view.ycam-_view.ycam1D;
  distz = _view.zcam-_view.zcam1D;
  edistx = _view.xeye-_view.xeye1D;
  edisty = _view.yeye-_view.yeye1D;
  edistz = _view.zeye-_view.zeye1D;
  ddistx = _view.dirx-_view.dirx1D;
  ddisty = _view.diry-_view.diry1D;
  ddistz = _view.dirz-_view.dirz1D;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx*froll;
    _view.ycam = _view.ycam - disty*froll;
    _view.zcam = _view.zcam - distz*froll;
    
    _view.xeye = _view.xeye - edistx*froll;
    _view.yeye = _view.yeye - edisty*froll;
    _view.zeye = _view.zeye - edistz*froll;

    _view.dirx = _view.dirx - ddistx*froll;
    _view.diry = _view.diry - ddisty*froll;
    _view.dirz = _view.dirz - ddistz*froll;

    display();
  }
}

//=============================================================================
/*
  Automatic movement to reach 2D saved camera position
  Current position is not saved
*/
//=============================================================================
void Data::rollto2Dws()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Compute distance
  distx = _view.xcam-_view.xcam2D;
  disty = _view.ycam-_view.ycam2D;
  distz = _view.zcam-_view.zcam2D;
  edistx = _view.xeye-_view.xeye2D;
  edisty = _view.yeye-_view.yeye2D;
  edistz = _view.zeye-_view.zeye2D;
  ddistx = _view.dirx-_view.dirx2D;
  ddisty = _view.diry-_view.diry2D;
  ddistz = _view.dirz-_view.dirz2D;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx*froll;
    _view.ycam = _view.ycam - disty*froll;
    _view.zcam = _view.zcam - distz*froll;
    
    _view.xeye = _view.xeye - edistx*froll;
    _view.yeye = _view.yeye - edisty*froll;
    _view.zeye = _view.zeye - edistz*froll;

    _view.dirx = _view.dirx - ddistx*froll;
    _view.diry = _view.diry - ddisty*froll;
    _view.dirz = _view.dirz - ddistz*froll;

    display();
  }
}

//=============================================================================
/*
  Automatic movement to reach 1D saved camera position
  Current position is not saved
*/
//=============================================================================
void Data::rollto1Dws()
{
  E_Int i;
  double distx, disty, distz;
  double edistx, edisty, edistz;
  double ddistx, ddisty, ddistz; 
  E_Int nroll = _pref.nroll;
  double froll = 1./((double)nroll);

  // Compute distance
  distx = _view.xcam - _view.xcam1D;
  disty = _view.ycam - _view.ycam1D;
  distz = _view.zcam - _view.zcam1D;
  edistx = _view.xeye - _view.xeye1D;
  edisty = _view.yeye - _view.yeye1D;
  edistz = _view.zeye - _view.zeye1D;
  ddistx = _view.dirx - _view.dirx1D;
  ddisty = _view.diry - _view.diry1D;
  ddistz = _view.dirz - _view.dirz1D;

  // Roll
  for (i = 0; i < nroll; i++)
  {
    _view.xcam = _view.xcam - distx*froll;
    _view.ycam = _view.ycam - disty*froll;
    _view.zcam = _view.zcam - distz*froll;
    
    _view.xeye = _view.xeye - edistx*froll;
    _view.yeye = _view.yeye - edisty*froll;
    _view.zeye = _view.zeye - edistz*froll;

    _view.dirx = _view.dirx - ddistx*froll;
    _view.diry = _view.diry - ddisty*froll;
    _view.dirz = _view.dirz - ddistz*froll;

    display();
  }
}

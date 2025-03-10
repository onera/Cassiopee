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

#include "Data.h"
#include <math.h>

// Global keyboard function
void gkeyboard(unsigned char key, int x, int y)
{
  Data* d = Data::getInstance();
  d->keyboard(key, x, y);
}

// Global keyboard up function
void gkeyboardup(int key, int x, int y)
{
  CPlotState* ptrState = Data::getInstance()->ptrState;
  if (ptrState->kkeysActivated == 0)
  {
    if (key == GLUT_KEY_UP) 
    {ptrState->keys[ptrState->kcursor]=5; ptrState->kcursor++;} 
    else if (key == GLUT_KEY_DOWN) 
    {ptrState->keys[ptrState->kcursor]=6; ptrState->kcursor++;}
    else if (key == GLUT_KEY_LEFT) 
    {ptrState->keys[ptrState->kcursor]=7; ptrState->kcursor++;}
    else if (key == GLUT_KEY_RIGHT) 
    {ptrState->keys[ptrState->kcursor]=8; ptrState->kcursor++;}
    return;
  }
}

// Global arrows function
void garrows(int key, int x, int y)
{
  Data* d = Data::getInstance();
  d->arrows(key, x, y);
}

//=============================================================================
// Keyboards calls
//=============================================================================
// Branche les differents appels pour le clavier
void Data::keyboard(unsigned char key, E_Int x, E_Int y)
{
  //printf("normal key: %d\n", key); fflush(stdout);
  E_Int nv;
  //E_Int stateHeader, stateInfo, stateMenu, stateBB;
  //double alpha = 0.04;
  //double dx = (_view.xeye-_view.xcam)*alpha;
  //double dy = (_view.yeye-_view.ycam)*alpha;
  //double dz = (_view.zeye-_view.zcam)*alpha;
  //double d = sqrt(dx*dx + dy*dy + dz*dz);
  //double dirx = _view.dirx;
  //double diry = _view.diry;
  //double dirz = _view.dirz;
  E_Int modif = glutGetModifiers();

  ptrState->render = 1;

  ptrState->keys[ptrState->kcursor] = key; ptrState->kcursor++;
  if (ptrState->kcursor > 127) ptrState->kcursor = 127;
  if (ptrState->kkeysActivated == 0) return; // no short cuts

  switch (key) {
  
  // -- Quit --
  case 'q':
  case 'Q':
  {
    glutHideWindow();
    freeGPUResources(-1, 0, _numberOfZones-1, 1);
    ptrState->freeGPURes = 1;
    _exit(0);
    break;
  }

  // -- Menu / display state information --
  // case 27: // esc key
  // {
  //   menu();
  //   break;
  // }
    
  // -- Move down --
  // case 'o':
  // {
  //   moveDown(alpha, dx, dy, dz, d, dirx, diry, dirz);
  //   break;
  // }
    
  // -- Move up --
  // case 'p':
  // {
  //   moveUp(alpha, dx, dy, dz, d, dirx, diry, dirz);
  //   break;
  // }
     
  // -- Fit view / fullscreen --
  case 'f':
  case 6: // on windows
  {
    if (modif == GLUT_ACTIVE_CTRL) 
    {
      if (ptrState->fullScreen == 0)
      { 
        glutFullScreen(); // pas de retour possible sous linux
        _view.wSav = _view.w; _view.hSav = _view.h;
        //glutReshapeWindow(1680, 1050);
        ptrState->fullScreen = 1;
      }
      else 
      {
        glutReshapeWindow(_view.wSav, _view.hSav);
        ptrState->fullScreen = 0;
      }
    }
    else { initCam(); farClipping(); }
    break;
  }
      
  // -- Mesh/Solid display --
  case '1':
  case 33: // !
  case 38: // 1
  {
    /*
    if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT))
    {
      if (ptrState->mode < 2) ptrState->mode++;
      else ptrState->mode = MESH;
    }
    else
    {
      if (ptrState->mode > 0) ptrState->mode--;
      else ptrState->mode = RENDER;
    }
    */

    if (ptrState->mode == SOLID) ptrState->mode = MESH;
    else if (ptrState->mode == MESH) ptrState->mode = SOLID;
    else ptrState->mode = MESH;
    break;
  }

  // -- Select/unselectAll --
  case 32: // space 
  {
    // toggle selectall
    E_Int toggle = 0;
    for (E_Int i = 0; i < _numberOfZones; i++)
      toggle = std::max(toggle, _zones[i]->selected);
    if (toggle == 1)
    {
      // unselect  all
      for (E_Int i = 0; i < _numberOfZones; i++) _zones[i]->selected = 0;
      ptrState->selectedZone = 0;
    }
    else
    {
      // select all
      for (E_Int i = 0; i < _numberOfZones; i++) _zones[i]->selected = 1;
    }
    break;
  }
    
  // -- Scalar mode --
  case '2':
  case 195: // 2
  case 233: // é
  {
    nv = _zones[0]->nfield;
    if (nv < 1) break;

    if (ptrState->mode != SCALARFIELD) 
    { ptrState->mode = SCALARFIELD; break; }
    
    if (modif == GLUT_ACTIVE_SHIFT)
    {
      ptrState->scalarField--;
      if (ptrState->scalarField < 0) ptrState->scalarField = nv-1;
    }
    else
    {
      ptrState->scalarField++;
      if (ptrState->scalarField >= nv) ptrState->scalarField = 0;
    }
    break;
  }
     
  // -- render mode --
  case '3':
  case 34: // 3
  {
    ptrState->mode = RENDER;
    break;
  }

  // -- I/J/K mode --
  case 'i':
  {
    changeIPlanePlus();
    break;
  }
  case 'I':
  {
    changeIPlaneMinus();
    break;
  }
  case 'j':
  {
    changeJPlanePlus();
    break;
  }
  case 'J':
  {
    changeJPlaneMinus();
    break;
  }
  case 'k':
  {
    changeKPlanePlus();
    break;
  }
  case 'K':
  {
    changeKPlaneMinus();
    break;
  }

  // -- Change selected zone in case of ambiguous selection --
  case '>':
  case '<':
  {
    changeAmbSelection();
    break;
  }

  // -- Change the dimension mode (3D - 2D - 1D) --
  case 'm':
  {
    switch (ptrState->dim)
    {
      case 3:
        ptrState->dim = 2;
        freeGPUResources(-1, 0, _numberOfZones-1, 0);
        roll3Dto2D();
        break;

      case 2:
        ptrState->dim = 3;
        freeGPUResources(-1, 0, _numberOfZones-1, 0);
        roll2Dto3D();
        break;
        
      case 1: // 1D mode est osolete dans cette version
        ptrState->dim = 3;
        roll1Dto3D();
        break;
    }
    break;
  }
  case 'M':
  {
    switch (ptrState->dim)
    {
      case 1:
        ptrState->dim = 2;
        roll1Dto2D();
        break;

      case 2:
        ptrState->dim = 3;
        freeGPUResources(-1, 0, _numberOfZones-1, 0);
        roll2Dto3D();
        break;
        
      case 3:
        ptrState->dim = 2;
        freeGPUResources(-1, 0, _numberOfZones-1, 0);
        roll3Dto2D();
        break;
    }
    break;
  }
      
  // Reload file
  case 'r':
  {
    PyEval_RestoreThread(_save);
    char com[1024];
    E_Int l = strlen(ptrState->file); char *p = ptrState->file;
    if (l > 5 && p[l-1] == 's' && p[l-2] == 'n' && p[l-3] == 'g' 
        && p[l-4] == 'c' && p[l-5] == '.')
      sprintf(com, "import Converter.PyTree; import CPlot.PyTree; kpl = Converter.PyTree.convertFile2PyTree('%s'); CPlot.PyTree.display(kpl)",
              ptrState->file);
    else
      sprintf(com, "import Converter; import CPlot; kpl = Converter.convertFile2Arrays('%s'); CPlot.display(kpl)",
              ptrState->file);
    PyRun_SimpleString(com);
    _save = PyEval_SaveThread(); 
    printTmpMessage("File reloaded.");
    break;
  }

  // Select active zone
  case 'z':
  {
    if (_pref.selectNextZone != NULL) _pref.selectNextZone->f(this);
    break;
  }
  case 'Z':
  {
    if (_pref.selectPreviousZone != NULL) _pref.selectPreviousZone->f(this);
    break;
  }

  // Look for active zone
  case 'l':
  {
    if (_pref.lookFor != NULL) _pref.lookFor->f(this);
    break;
  }

  // Zone deactivation/reactivation
  case 'a':
  {
    if (ptrState->selectedZone != 0)
    {
      _zones[ptrState->selectedZone-1]->active = 0;
      _zones[ptrState->selectedZone-1]->selected = 0;
      /*
      if (ptrState->deactivatedZones == NULL)
      {
        struct chain_int* ci;
        ci = (struct chain_int*)malloc(sizeof(struct chain_int));
        ci->value = ptrState->selectedZone;
        ci->next = NULL;
        ptrState->deactivatedZones = ci;
      }
      else
      {
        struct chain_int* ci = ptrState->deactivatedZones;
        while (ci->next != NULL) ci = ci->next;
        ci->next = (struct chain_int*)malloc(sizeof(struct chain_int));
        ci = ci->next;
        ci->value = ptrState->selectedZone;
        ci->next = NULL;
      } */
      ptrState->insertDeactivatedZones(ptrState->selectedZone);
      //ptrState->printDeactivatedZones();

      ptrState->selectedZone = ptrState->selectedZone+1;
      if (ptrState->selectedZone == _numberOfZones+1) 
        ptrState->selectedZone = 0;
      else
      {
        while (_zones[ptrState->selectedZone-1]->active == 0)
        {
          ptrState->selectedZone = ptrState->selectedZone+1;
          if (ptrState->selectedZone == _numberOfZones+1) 
          {
            ptrState->selectedZone = 0;
            break;
          }
        }
        if (ptrState->selectedZone != 0)
          _zones[ptrState->selectedZone-1]->selected = 1;
      }
    }
    break;
  }
  case 'A':
  {
    if (ptrState->deactivatedZones != NULL)
    {
      struct chain_int* ci = ptrState->deactivatedZones;
      struct chain_int* cip = NULL;
      while (ci->next != NULL)
      {
        cip = ci;
        ci = ci->next;
      }
      if (ptrState->selectedZone != 0) _zones[ptrState->selectedZone-1]->selected = 0;
      _zones[ci->value-1]->active = 1;
      //ptrState->selectedZone = ci->value;
      if (ptrState->selectedZone != 0) _zones[ptrState->selectedZone-1]->selected = 1;
      if (cip == NULL) ptrState->deactivatedZones = NULL;
      else cip->next = NULL;
      free(ci);
    }
    //ptrState->printDeactivatedZones();
    break;
  }

  // Image dump
  // case 'y':
  // {
  //   stateHeader = ptrState->header;
  //   stateInfo = ptrState->info;
  //   stateMenu = ptrState->menu;
  //   stateBB = ptrState->bb;
  //   ptrState->header = 0;
  //   ptrState->info = 0;
  //   ptrState->menu = 0;
  //   ptrState->bb = 0;
  //   display();
  //   dumpWindow();
  //   ptrState->header = stateHeader;
  //   ptrState->info = stateInfo;
  //   ptrState->menu = stateMenu;
  //   ptrState->bb = stateBB;
  //   printTmpMessage("Image dumped to file.");
  //   break;
  // }

  // Change render appearance
  case 'c':
  {
    changeAppearance();
    break;
  }

  default:
  break;

  }
}

//=============================================================================
// Branche les appels suivants les fleches
//=============================================================================
void Data::arrows(unsigned char key, E_Int x, E_Int y)
{
  double alpha = 0.07853981633974483;
  double dx = (_view.xeye - _view.xcam)*alpha;
  double dy = (_view.yeye - _view.ycam)*alpha;
  double dz = (_view.zeye - _view.zcam)*alpha;
  double d = sqrt(dx*dx + dy*dy + dz*dz);
  E_Int modif = glutGetModifiers();
  double dirx = _view.dirx;
  double diry = _view.diry;
  double dirz = _view.dirz;

  ptrState->render = 1;
  //printf("%d\n", key);
  
  if (ptrState->kkeysActivated == 0)
  {
    //printf("special key %d\n", key);
    if (key == GLUT_KEY_UP) 
    {ptrState->keys[ptrState->kcursor]=1; ptrState->kcursor++;} 
    else if (key == GLUT_KEY_DOWN) 
    {ptrState->keys[ptrState->kcursor]=2; ptrState->kcursor++;}
    else if (key == GLUT_KEY_LEFT) 
    {ptrState->keys[ptrState->kcursor]=3; ptrState->kcursor++;}
    else if (key == GLUT_KEY_RIGHT) 
    {ptrState->keys[ptrState->kcursor]=4; ptrState->kcursor++;}
    else 
    {ptrState->keys[ptrState->kcursor]=key; ptrState->kcursor++;}
    return;
  }

  switch (key)
  {
    case GLUT_KEY_UP:
    {
      if (modif == GLUT_ACTIVE_SHIFT)
        strafeUp(alpha, dx, dy, dz, d, dirx, diry, dirz);
      else if (modif == GLUT_ACTIVE_CTRL)
      {
        _view.xcam += dx;
        _view.ycam += dy;
        _view.zcam += dz;
        adaptiveClipping(d);
      }
      else
        moveUp(alpha, dx, dy, dz, d, dirx, diry, dirz);
      break;
    }
      
    case GLUT_KEY_DOWN:
    {
      if (modif == GLUT_ACTIVE_SHIFT)
        strafeDown(alpha, dx, dy, dz, d, dirx, diry, dirz);
      else if (modif == GLUT_ACTIVE_CTRL)
      {
        _view.xcam -= dx;
        _view.ycam -= dy;
        _view.zcam -= dz;
        adaptiveClipping(d);
      }
      else
        moveDown(alpha, dx, dy, dz, d, dirx, diry, dirz);
      break;
    }
      
    case GLUT_KEY_LEFT:
    {
      if (modif == GLUT_ACTIVE_SHIFT)
        strafeRight(alpha, dx, dy, dz, d, dirx, diry, dirz);
      else if (modif == GLUT_ACTIVE_CTRL)
        tiltLeft(alpha, dx, dy, dz, d, dirx, diry, dirz);
      else
        moveRight(alpha, dx, dy, dz, d, dirx, diry, dirz);
      break;
    }
      
    case GLUT_KEY_RIGHT:
    {
      if (modif == GLUT_ACTIVE_SHIFT)
        strafeLeft(alpha, dx, dy, dz, d, dirx, diry, dirz);
      else if (modif == GLUT_ACTIVE_CTRL)
        tiltRight(alpha, dx, dy, dz, d, dirx, diry, dirz);
      else
        moveLeft(alpha, dx, dy, dz, d, dirx, diry, dirz);
      break;
    }
  }
  //printf("camera position %f %f %f\n",_view.xcam,_view.ycam,_view.zcam);
}

//=============================================================================
// Move down
//=============================================================================
void Data::moveDown(double alpha, double dx, double dy, double dz, double d,
                    double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    double P0ex, P0ey, P0ez, ox, oy, oz, P1x, P1y, P1z;
    double P1ex, P1ey, P1ez, nv, nd, f, nP1e, nP0e, k;
    P1x = _view.xcam - dirx*d;
    P1y = _view.ycam - diry*d;
    P1z = _view.zcam - dirz*d;
    P1ex = P1x - _view.xeye;
    P1ey = P1y - _view.yeye;
    P1ez = P1z - _view.zeye;

    P0ex = _view.xcam - _view.xeye;
    P0ey = _view.ycam - _view.yeye;
    P0ez = _view.zcam - _view.zeye;

    ox = P0ey*dirz - P0ez*diry;
    oy = P0ez*dirx - P0ex*dirz;
    oz = P0ex*diry - P0ey*dirx;
    nv = dirx*dirx+diry*diry+dirz*dirz;
    
    dx = P1ey*oz - P1ez*oy;
    dy = P1ez*ox - P1ex*oz;
    dz = P1ex*oy - P1ey*ox;

    nd = dx*dx + dy*dy + dz*dz;
    if (nd > 1.e-24) f = sqrt(nv/nd);
    else f = 1.;
    dx = - dx * f;
    dy = - dy * f;
    dz = - dz * f;

    _view.dirx = dx;
    _view.diry = dy;
    _view.dirz = dz;

    nP1e = P1ex*P1ex+P1ey*P1ey+P1ez*P1ez;
    nP0e = P0ex*P0ex+P0ey*P0ey+P0ez*P0ez;
    if (nP1e > 1.e-24) k = sqrt(nP0e / nP1e);
    else k = 0.;    
    _view.xcam = _view.xeye + k*P1ex;
    _view.ycam = _view.yeye + k*P1ey;
    _view.zcam = _view.zeye + k*P1ez;

  }
  else
  {
    // 2D mode
    switch (ptrState->var2D)
    {
      case 0:
        _view.ycam = _view.ycam - d;
        _view.yeye = _view.yeye - d;
      break;

      case 1:
        _view.zcam = _view.zcam - d;
        _view.zeye = _view.zeye - d;
        break;

      case 2:
        _view.zcam = _view.zcam - d;
        _view.zeye = _view.zeye - d;
        break;
    }
  }
}

//=============================================================================
// Strafe down
//=============================================================================
void Data::strafeDown(double alpha, double dx, double dy, double dz, double d,
                      double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xcam = _view.xcam - dirx*d;
    _view.ycam = _view.ycam - diry*d;
    _view.zcam = _view.zcam - dirz*d;
    _view.xeye = _view.xeye - dirx*d;
    _view.yeye = _view.yeye - diry*d;
    _view.zeye = _view.zeye - dirz*d;
  }
  else
  {
    // 2D mode
    switch (ptrState->var2D)
    {
      case 0:
        _view.ycam = _view.ycam - d;
        _view.yeye = _view.yeye - d;
      break;

      case 1:
        _view.zcam = _view.zcam - d;
        _view.zeye = _view.zeye - d;
        break;

      case 2:
        _view.zcam = _view.zcam - d;
        _view.zeye = _view.zeye - d;
        break;
    }
  }
}

//=============================================================================
// Move up
//=============================================================================
void Data::moveUp(double alpha, double dx, double dy, double dz, double d,
                  double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    double P0ex, P0ey, P0ez, ox, oy, oz, P1x, P1y, P1z;
    double P1ex, P1ey, P1ez, nv, nd, f, nP1e, nP0e, k;
    P1x = _view.xcam + dirx*d;
    P1y = _view.ycam + diry*d;
    P1z = _view.zcam + dirz*d;
    P1ex = P1x - _view.xeye;
    P1ey = P1y - _view.yeye;
    P1ez = P1z - _view.zeye;
 
    P0ex = _view.xcam - _view.xeye;
    P0ey = _view.ycam - _view.yeye;
    P0ez = _view.zcam - _view.zeye;

    ox = P0ey*dirz - P0ez*diry;
    oy = P0ez*dirx - P0ex*dirz;
    oz = P0ex*diry - P0ey*dirx;
    nv = dirx*dirx+diry*diry+dirz*dirz;
    
    dx = P1ey*oz - P1ez*oy;
    dy = P1ez*ox - P1ex*oz;
    dz = P1ex*oy - P1ey*ox;

    nd = dx*dx + dy*dy + dz*dz;
    if (nd > 1.e-24) f = sqrt(nv/nd);
    else f = 1.;
    dx = - dx * f;
    dy = - dy * f;
    dz = - dz * f;

    _view.dirx = dx;
    _view.diry = dy;
    _view.dirz = dz;

    nP1e = P1ex*P1ex+P1ey*P1ey+P1ez*P1ez;
    nP0e = P0ex*P0ex+P0ey*P0ey+P0ez*P0ez;
    if (nP1e > 1.e-24) k = sqrt(nP0e / nP1e);
    else k = 0.;    
    _view.xcam = _view.xeye + k*P1ex;
    _view.ycam = _view.yeye + k*P1ey;
    _view.zcam = _view.zeye + k*P1ez;

  }
  else
  {
    // 2D mode
    switch (ptrState->var2D)
    {
      case 0:
        _view.ycam = _view.ycam + d;
        _view.yeye = _view.yeye + d;
        break;
        
      case 1:
        _view.zcam = _view.zcam + d;
        _view.zeye = _view.zeye + d;
        break;

      case 2:
        _view.zcam = _view.zcam + d;
        _view.zeye = _view.zeye + d;
        break;
    }
  }
}

//=============================================================================
// Strafe up
//=============================================================================
void Data::strafeUp(double alpha, double dx, double dy, double dz, double d,
                    double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xcam = _view.xcam + dirx*d;
    _view.ycam = _view.ycam + diry*d;
    _view.zcam = _view.zcam + dirz*d;
    _view.xeye = _view.xeye + dirx*d;
    _view.yeye = _view.yeye + diry*d;
    _view.zeye = _view.zeye + dirz*d;
  }
  else
  {
    // 2D mode
    switch (ptrState->var2D)
    {
      case 0:
        _view.ycam = _view.ycam + d;
        _view.yeye = _view.yeye + d;
        break;
        
      case 1:
        _view.zcam = _view.zcam + d;
        _view.zeye = _view.zeye + d;
        break;

      case 2:
        _view.zcam = _view.zcam + d;
        _view.zeye = _view.zeye + d;
        break;
    }
  }
}

//=============================================================================
// Move right
//=============================================================================
void Data::moveRight(double alpha, double dx, double dy, double dz, double d,
                     double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    /*
    _view.xcam = _view.xcam + dy*dirz - dz*diry;
    _view.ycam = _view.ycam - dx*dirz + dz*dirx;
    _view.zcam = _view.zcam + dx*diry - dy*dirx;
    */
    double P1x,P1y,P1z,P1ex,P1ey,P1ez,P0ex,P0ey,P0ez,nP1e,nP0e,k;
    P1x = _view.xcam + dy*dirz - dz*diry;
    P1y = _view.ycam - dx*dirz + dz*dirx;
    P1z = _view.zcam + dx*diry - dy*dirx;
    P1ex = P1x - _view.xeye;
    P1ey = P1y - _view.yeye;
    P1ez = P1z - _view.zeye;
    P0ex = _view.xcam - _view.xeye;
    P0ey = _view.ycam - _view.yeye;
    P0ez = _view.zcam - _view.zeye;
    nP1e = P1ex*P1ex+P1ey*P1ey+P1ez*P1ez;
    nP0e = P0ex*P0ex+P0ey*P0ey+P0ez*P0ez;
    if (nP1e > 1.e-24) k = sqrt(nP0e / nP1e);
    else k = 0.;
    _view.xcam = _view.xeye + k*P1ex;
    _view.ycam = _view.yeye + k*P1ey;
    _view.zcam = _view.zeye + k*P1ez;
  }
  else if (ptrState->dim == 2)
  {
    switch ( ptrState->var2D )
    {
      case 0:
        _view.xcam = _view.xcam + d;
        _view.xeye = _view.xeye + d;
        break;
      case 1:
        _view.xcam = _view.xcam + d;
        _view.xeye = _view.xeye + d;
        break;
      case 2:
        _view.ycam = _view.ycam + d;
        _view.yeye = _view.yeye + d;
        break;
    }
  }
  else
  {
    switch ( ptrState->var1D )
    {
      case 0:
        _view.xcam = _view.xcam + d;
        _view.xeye = _view.xeye + d;
        break;
      case 1:
        _view.ycam = _view.ycam + d;
        _view.yeye = _view.yeye + d;
        break;
      case 2:
        _view.zcam = _view.zcam + d;
        _view.zeye = _view.zeye + d;
        break;
    }
  }
}

//=============================================================================
// Strafe right
//=============================================================================
void Data::strafeRight(double alpha, double dx, double dy, double dz, double d,
                       double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xcam = _view.xcam + dy*dirz - dz*diry;
    _view.ycam = _view.ycam - dx*dirz + dz*dirx;
    _view.zcam = _view.zcam + dx*diry - dy*dirx;
    _view.xeye = _view.xeye + dy*dirz - dz*diry;
    _view.yeye = _view.yeye - dx*dirz + dz*dirx;
    _view.zeye = _view.zeye + dx*diry - dy*dirx;
  }
  else if (ptrState->dim == 2)
  {
    switch ( ptrState->var2D )
    {
      case 0:
        _view.xcam = _view.xcam + d;
        _view.xeye = _view.xeye + d;
        break;
      case 1:
        _view.xcam = _view.xcam + d;
        _view.xeye = _view.xeye + d;
        break;
      case 2:
        _view.ycam = _view.ycam - d;
        _view.yeye = _view.yeye - d;
        break;
    }
  }
  else
  {
    switch ( ptrState->var1D )
    {
      case 0:
        _view.xcam = _view.xcam + d;
        _view.xeye = _view.xeye + d;
        break;
      case 1:
        _view.ycam = _view.ycam + d;
        _view.yeye = _view.yeye + d;
        break;
      case 2:
        _view.zcam = _view.zcam + d;
        _view.zeye = _view.zeye + d;
        break;
    }
  }
}

//=============================================================================
// Tilt left
//=============================================================================
void Data::tiltLeft(double alpha, double dx, double dy, double dz, double d,
                    double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    double z1, z2, z3, d1, d2, d3, n;
    z1 = dy*dirz - dz*diry;
    z2 = dz*dirx - dx*dirz;
    z3 = dx*diry - dy*dirx;
    n = 1./sqrt(z1*z1 + z2*z2 + z3*z3);
    z1 = z1*n; z2 = z2*n; z3 = z3*n;
    d1 = dirx + alpha*z1;
    d2 = diry + alpha*z2;
    d3 = dirz + alpha*z3;
    z1 = 1./sqrt(d1*d1 + d2*d2 + d3*d3);
    _view.dirx =  d1 * z1;
    _view.diry =  d2 * z1;
    _view.dirz =  d3 * z1;
  }
}

//=============================================================================
// Tilt right
//=============================================================================
void Data::tiltRight(double alpha, double dx, double dy, double dz, double d,
                     double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
     double z1, z2, z3, d1, d2, d3, n;
    z1 = dy*dirz - dz*diry;
    z2 = dz*dirx - dx*dirz;
    z3 = dx*diry - dy*dirx;
    n = 1./sqrt(z1*z1 + z2*z2 + z3*z3);
    z1 = z1*n; z2 = z2*n; z3 = z3*n;
    d1 = dirx - alpha*z1;
    d2 = diry - alpha*z2;
    d3 = dirz - alpha*z3;
    z1 = 1./sqrt(d1*d1 + d2*d2 + d3*d3);
    _view.dirx =  d1 * z1;
    _view.diry =  d2 * z1;
    _view.dirz =  d3 * z1;
  }
}

//=============================================================================
// Rotate head right
//=============================================================================
void Data::rotateHeadRight(double alpha, double dx, double dy, double dz, 
                           double d,
                           double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xeye = _view.xeye - dy*dirz + dz*diry;
    _view.yeye = _view.yeye + dx*dirz - dz*dirx;
    _view.zeye = _view.zeye - dx*diry + dy*dirx; 
  }
}

//=============================================================================
// Rotate head left
//=============================================================================
void Data::rotateHeadLeft(double alpha, double dx, double dy, double dz, 
                          double d,
                          double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xeye = _view.xeye + dy*dirz - dz*diry;
    _view.yeye = _view.yeye - dx*dirz + dz*dirx;
    _view.zeye = _view.zeye + dx*diry - dy*dirx; 
  }
}

//=============================================================================
// Rotate head up
//=============================================================================
void Data::rotateHeadUp(double alpha, double dx, double dy, double dz, 
                        double d,
                        double dirx, double diry, double dirz)
{ 
  if (ptrState->dim == 3)
  {
    _view.xeye = _view.xeye + d*dirx;
    _view.yeye = _view.yeye + d*diry;
    _view.zeye = _view.zeye + d*dirz; 
  }
}

//=============================================================================
// Rotate down up
//=============================================================================
void Data::rotateHeadDown(double alpha, double dx, double dy, double dz, 
                          double d,
                          double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xeye = _view.xeye - d*dirx;
    _view.yeye = _view.yeye - d*diry;
    _view.zeye = _view.zeye - d*dirz; 
  }
}

//=============================================================================
// Move left
//=============================================================================
void Data::moveLeft(double alpha, double dx, double dy, double dz, double d,
                    double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    /*
    _view.xcam = _view.xcam - dy*dirz + dz*diry;
    _view.ycam = _view.ycam + dx*dirz - dz*dirx;
    _view.zcam = _view.zcam - dx*diry + dy*dirx;
    */
    double P1x,P1y,P1z,P1ex,P1ey,P1ez,P0ex,P0ey,P0ez,nP1e,nP0e,k;
    P1x = _view.xcam - dy*dirz + dz*diry;
    P1y = _view.ycam + dx*dirz - dz*dirx;
    P1z = _view.zcam - dx*diry + dy*dirx;
    P1ex = P1x - _view.xeye;
    P1ey = P1y - _view.yeye;
    P1ez = P1z - _view.zeye;
    P0ex = _view.xcam - _view.xeye;
    P0ey = _view.ycam - _view.yeye;
    P0ez = _view.zcam - _view.zeye;
    nP1e = P1ex*P1ex+P1ey*P1ey+P1ez*P1ez;
    nP0e = P0ex*P0ex+P0ey*P0ey+P0ez*P0ez;
    if (nP1e > 1.e-24) k = sqrt(nP0e / nP1e);
    else k = 0.;
    _view.xcam = _view.xeye + k*P1ex;
    _view.ycam = _view.yeye + k*P1ey;
    _view.zcam = _view.zeye + k*P1ez;
  }
  else if (ptrState->dim == 2)
  {
    switch (ptrState->var2D)
    {
      case 0:
        _view.xcam = _view.xcam - d;
        _view.xeye = _view.xeye - d;
        break;
      case 1:
        _view.xcam = _view.xcam - d;
        _view.xeye = _view.xeye - d;
        break;
      case 2:
        _view.ycam = _view.ycam - d;
        _view.yeye = _view.yeye - d;
        break;
    }
  }
  else
  {
    switch (ptrState->var1D)
    {
      case 0:
        _view.xcam = _view.xcam - d;
        _view.xeye = _view.xeye - d;
        break;
      case 1:
        _view.ycam = _view.ycam - d;
        _view.yeye = _view.yeye - d;
        break;
      case 2:
        _view.zcam = _view.zcam - d;
        _view.zeye = _view.zeye - d;
        break;
    }
  }
}

//=============================================================================
// Strafe left
//=============================================================================
void Data::strafeLeft(double alpha, double dx, double dy, double dz, double d,
                      double dirx, double diry, double dirz)
{
  if (ptrState->dim == 3)
  {
    _view.xcam = _view.xcam - dy*dirz + dz*diry;
    _view.ycam = _view.ycam + dx*dirz -dz*dirx;
    _view.zcam = _view.zcam - dx*diry + dy*dirx; 
    _view.xeye = _view.xeye - dy*dirz + dz*diry;
    _view.yeye = _view.yeye + dx*dirz -dz*dirx;
    _view.zeye = _view.zeye - dx*diry + dy*dirx; 
  }
  else if (ptrState->dim == 2)
  {
    switch (ptrState->var2D)
    {
      case 0:
        _view.xcam = _view.xcam - d;
        _view.xeye = _view.xeye - d;
        break;
      case 1:
        _view.xcam = _view.xcam - d;
        _view.xeye = _view.xeye - d;
        break;
      case 2:
        _view.ycam = _view.ycam + d;
        _view.yeye = _view.yeye + d;
        break;
    }
  }
  else
  {
    switch (ptrState->var1D)
    {
      case 0:
        _view.xcam = _view.xcam - d;
        _view.xeye = _view.xeye - d;
        break;
      case 1:
        _view.ycam = _view.ycam - d;
        _view.yeye = _view.yeye - d;
        break;
      case 2:
        _view.zcam = _view.zcam - d;
        _view.zeye = _view.zeye - d;
        break;
    }
  }
}

//=============================================================================
// Change i plane plus
//=============================================================================
void Data::changeIPlanePlus()
{
  E_Int modif = glutGetModifiers();
  E_Int* changed = new E_Int [_numberOfStructZones];
  E_Int nchanged = 0;

  for (E_Int nz = 0; nz < _numberOfStructZones; nz++)
  {
    StructZone* z = _szones[nz];
    if (z->selected == 1)
    {
      z->activePlane = 0;
      z->iPlane++;
      if (z->iPlane > z->ni-1) z->iPlane = -1;
      if (modif == GLUT_ACTIVE_CTRL) z->iPlane = -1;
      if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) z->iPlane = -2; // no plane display
      changed[nchanged] = nz;
      nchanged++;
      z->destroyGPUIsoField();
    }
  }
  updateGPUResources( -1, nchanged, 0, (void*)changed );
}

//=============================================================================
// Change j plane plus
//=============================================================================
void Data::changeJPlanePlus()
{
  E_Int modif = glutGetModifiers();
  E_Int* changed = new E_Int [_numberOfStructZones];
  E_Int nchanged = 0;

  for (E_Int nz = 0; nz < _numberOfStructZones; nz++)
  {
    StructZone* z = _szones[nz];
    if (z->selected == 1)
    {
      z->activePlane = 1;
      z->jPlane++;
      if (z->jPlane > z->nj-1) z->jPlane = -1;
      if (modif == GLUT_ACTIVE_CTRL) z->jPlane = -1;
      if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) z->jPlane = -2; // no plane display
      changed[nchanged] = nz;
      nchanged++;
      z->destroyGPUIsoField();
    }
  }
  updateGPUResources( -1, nchanged, 0, (void*)changed );
}

//=============================================================================
// Change k plane plus
//=============================================================================
void Data::changeKPlanePlus()
{
  E_Int modif = glutGetModifiers();
  E_Int* changed = new E_Int [_numberOfStructZones];
  E_Int nchanged = 0;

  for (E_Int nz = 0; nz < _numberOfStructZones; nz++)
  {
    StructZone* z = _szones[nz];
    if (z->selected == 1)
    {
      z->activePlane = 2;
      z->kPlane++;
      if (z->kPlane > z->nk-1) z->kPlane = -1;
      if (modif == GLUT_ACTIVE_CTRL) z->kPlane = -1; // both plane display
      if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) z->kPlane = -2; // no plane display
      changed[nchanged] = nz;
      nchanged++;
      z->destroyGPUIsoField();
    }
  }
  updateGPUResources(-1, nchanged, 0, (void*)changed);
}

//=============================================================================
// Change i plane minus
//=============================================================================
void Data::changeIPlaneMinus()
{
  E_Int modif = glutGetModifiers();
  E_Int* changed = new E_Int [_numberOfStructZones];
  E_Int nchanged = 0;

  for (E_Int nz = 0; nz < _numberOfStructZones; nz++)
  {
    StructZone* z = _szones[nz];
    if (z->selected == 1)
    {
      z->activePlane = 0;
      z->iPlane--;
      if (z->iPlane < -1) z->iPlane = z->ni-1;
      if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) z->iPlane = -2; // no plane display
      changed[nchanged] = nz;
      nchanged++;
      z->destroyGPUIsoField();
    }
  }
  updateGPUResources(-1, nchanged, 0, (void*)changed);
}

//=============================================================================
// Change j plane minus
//=============================================================================
void Data::changeJPlaneMinus()
{
  E_Int modif = glutGetModifiers();
  E_Int* changed = new E_Int [_numberOfStructZones];
  E_Int nchanged = 0;

  for (E_Int nz = 0; nz < _numberOfStructZones; nz++)
  {
    StructZone* z = _szones[nz];
    if (z->selected == 1)
    {
      z->activePlane = 1;
      z->jPlane--;
      if (z->jPlane < -1) z->jPlane = z->nj-1;
      if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) z->jPlane = -2; // no plane display
      changed[nchanged] = nz;
      nchanged++;
      z->destroyGPUIsoField();
    }
  }
  updateGPUResources(-1, nchanged, 0, (void*)changed);
}

//=============================================================================
// Change k plane minus
//=============================================================================
void Data::changeKPlaneMinus()
{
  E_Int modif = glutGetModifiers();
  E_Int* changed = new E_Int [_numberOfStructZones];
  E_Int nchanged = 0;

  for (E_Int nz = 0; nz < _numberOfStructZones; nz++)
  {
    StructZone* z = _szones[nz];
    if (z->selected == 1)
    {
      z->activePlane = 1;
      z->kPlane--;
      if (z->kPlane < -1) z->kPlane = z->nk-1;
      if (modif == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) z->kPlane = -2; // no plane display
      changed[nchanged] = nz;
      nchanged++;
      z->destroyGPUIsoField();
    }
  }
  updateGPUResources(-1, nchanged, 0, (void*)changed);
}

//=============================================================================
// Change secondary variable plus
//=============================================================================
void Data::changeSecondaryVariablePlus()
{
  if (ptrState->dim == 2)
  {
    ptrState->var2D++;
    if (ptrState->var2D > 2) ptrState->var2D = 0;
    init2DCam();
    rollto2Dws();
  }
  else
  {
    ptrState->var1D++;
    if (ptrState->var1D > 3) ptrState->var1D = 0;
    init1DCam();
    rollto1Dws();
  }
}

//=============================================================================
// Change secondary variable plus
//=============================================================================
void Data::changeSecondaryVariableMinus()
{
  if (ptrState->dim == 2)
  {
    ptrState->var2D--;
    if (ptrState->var2D < 0) ptrState->var2D = 2;
    init2DCam();
    rollto2Dws();
  }
  else
  {
    ptrState->var1D--;
    if (ptrState->var1D < 0) ptrState->var1D = 3;
    init1DCam();
    rollto1Dws();
  }
}

//=============================================================================
// Change blanking function
//=============================================================================
void Data::changeBlankingFunction()
{
  if (_pref.blanking == NULL)
  {
    _pref.blanking = _plugins.blanking;
    if (_pref.blanking == NULL)
    {
      printf("no blanking avaiblable\n");
      printTmpMessage("No blanking available.");
    }
    else
    {
      findBlankedZones();
      printTmpMessage(_pref.blanking->functionName);
    }
  }
  else
  {
    _pref.blanking = _pref.blanking->next;
    findBlankedZones();
    
    if (_pref.blanking == NULL)
      printTmpMessage("No blanking set.");
    else
      printTmpMessage(_pref.blanking->functionName);
  }
}

//=============================================================================
// Change the render appearance 
//=============================================================================
void Data::changeAppearance()
{
  // // In iso solid mode, change the colormap or the light
  // if (ptrState->mode == SCALARFIELD)
  // {
  //   if (ptrState->isoLight == 0)
  //   {
  //     ptrState->isoLight = 1;
  //     //printTmpMessage("Activating light for iso.");
  //   }
  //   else
  //   {
  //     if (_pref.colorMap->next == NULL)
  //       _pref.colorMap = _plugins.colorMap;
  //     else
  //       _pref.colorMap = _pref.colorMap->next;
  //     //printTmpMessage(_pref.colorMap->functionName);
  //     ptrState->isoLight = 0;
  //   }
  // }

  // In solid mode 
  if (ptrState->mode == SOLID)
  {
    ptrState->solidStyle = ptrState->solidStyle+1;
    if (ptrState->solidStyle == 2) ptrState->solidStyle = 3; // bypass white
    if (ptrState->solidStyle > 4) ptrState->solidStyle = 0;
  }

  // In mesh mode
  if (ptrState->mode == MESH)
  {
    ptrState->meshStyle = ptrState->meshStyle+1;
    if (ptrState->meshStyle > 4) ptrState->meshStyle = 0;
  }
}

//=============================================================================
// Change ambiguous selection 
//=============================================================================
void Data::changeAmbSelection()
{
  std::vector<E_Int>& v = ptrState->ambSelections;
  E_Int zone, ind, inde, ncon;
  E_Int vsize = v.size()/4;
  if (vsize <= 1) return;
  Zone* z = NULL;

  // find next active zone in ambiguous selections
  //printf("already: %d\n", ptrState->ambSelSet);
  E_Int roll = 0;
  while (roll < vsize)
  {
    ptrState->ambSelSet += 1;
    if (ptrState->ambSelSet >= vsize) ptrState->ambSelSet = 0;
    zone = v[ptrState->ambSelSet*4];
    z = _zones[zone];
    if (z->active == 1) break;
    roll++;
  }
  if (z->active == 0) return; // no available active zone
  //printf("new: %d\n", ptrState->ambSelSet);
  
  ind = v[ptrState->ambSelSet*4+1];
  inde = v[ptrState->ambSelSet*4+2];
  ncon = v[ptrState->ambSelSet*4+3];
  
  for (E_Int i = 0; i < _numberOfZones; i++)
  _zones[i]->previouslySelected = _zones[i]->selected;

  ptrState->selectedZone = zone+1;

  for (E_Int i = 0; i < _numberOfZones; i++) _zones[i]->selected = 0;
  z->selected = 1;

  if (zone < _numberOfStructZones)
  {
    StructZone* zz = (StructZone*)z;
    E_Int ni = zz->ni; 
    E_Int nj = zz->nj;
    E_Int k = ind / (ni*nj);
    E_Int j = (ind - k*ni*nj)/ni;
    E_Int i = ind - k*ni*nj - j*ni;
    ptrState->activePointI = i+1;
    ptrState->activePointJ = j+1;
    ptrState->activePointK = k+1;
  }
  else
  {
    ptrState->activePointI = ind; // indice du noeud le plus proche
    ptrState->activePointJ = inde; // indice de l'element contenant P
    ptrState->activePointL = ncon; // connectivite contenant l'element
    UnstructZone* zz = (UnstructZone*)z;
    if (zz->eltType[0] != 10) // autre que NGON
    {
      E_Int* c = zz->connect[ncon];
      E_Int size = zz->eltSize[ncon];
      E_Int ne = zz->nec[ncon];
      E_Int v = 0;
      E_Int prev = 0;
      for (E_Int nc = 0; nc < ncon; nc++) prev += zz->nec[nc];
      for (E_Int nv = 0; nv < size; nv++)
      {
        if (c[inde-prev+nv*ne] == ind+1) break;
      }
      ptrState->activePointK = -v-1;
    }
    //else ptrState->activePointK = findFace(posX, posY, posZ, indE, zz, dist);
  }
  for (E_Int n = 0; n < z->nfield; n++)
  {
    double* f = z->f[n];
    ptrState->activePointF[n] = f[ind];
  }

}

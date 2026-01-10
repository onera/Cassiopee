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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//=============================================================================
/*
  Display une zone en iso solides pour les champs scalaires.
  IN: zonep: pointeur sur la zone a afficher
  IN: zone: le no de la zone dans la liste globale des zones
  IN: nofield: le no du champ
*/
//=============================================================================
void DataDL::renderSIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield)
{
  //E_Int i, j, k, n1, n2;
  float offb;
  double blend;
  
  // Grid dimensions
  E_Int ni = zonep->ni;
  E_Int nj = zonep->nj;
  E_Int nk = zonep->nk;
  
  ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  // Blending
  blend = 1.;
#include "selection2.h"
  bool is1D = false;
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1) is1D = true;
  

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend

  if (is1D)
  {
    if (curr != 0) _shaders[curr]->setUniform("lightOn", (int)0); // impose isoLight off on 1D meshes
  }
#endif

  glCallList(zoneImpl->_DLiso);
      
#ifdef __SHADERS__
  if (is1D)
  {
    if (ptrState->isoLight == 1 && ptrState->dim == 3)
    {
      if (curr != 0) _shaders[curr]->setUniform("lightOn", (int)1); // put back the isoLight value found in the CPlot state
    }
  }
#endif
}

//=============================================================================
/*
  Display une zone en iso solides pour les champs vectoriels.
  IN: zonep: pointeur sur la zone a afficher
  IN: zone: le no de la zone dans la liste globale des zones
  IN: nofield1, nofield2, nofield3: les no des champs
*/
//=============================================================================
void DataDL::renderSIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield1,
                                 E_Int nfield2, E_Int nofield3)
{  
  E_Int i, j, k, n1, n2;
  float offb;
  double blend;
  E_Int ret1, ret2;

  // Grid dimensions
  E_Int ni = zonep->ni;
  E_Int nj = zonep->nj;
  E_Int nk = zonep->nk;
  if (ptrState->dim == 2) nk = 1;
  E_Int nij = ni*nj;

  ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  // Blending
  blend = 1.;
#include "selection2.h"
  bool is1D = false;
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1) is1D = true;

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

  double* x = zonep->x; double* y = zonep->y; double* z = zonep->z;

  glCallList(zoneImpl->_DLiso);
      
  // Pour les lignes
  if (is1D)
  {
    glBegin(GL_LINES);
    E_Int nie, nje, nke;
    nie = ni; nje = nj; nke = nk;
    if (ni*nj == 1) nke = nke-1;
    if (ni*nk == 1) nje = nje-1;
    if (nj*nk == 1) nie = nie-1;
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            glColor4f(0, 0, 0+offb, blend);
            glVertex3d(x[n1], y[n1], z[n1]);
            glColor4f(0, 0, 0+offb, blend);
            glVertex3d(x[n2], y[n2], z[n2]);
          }
    }
    else
    {
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              glColor4f(0, 0, 0+offb, blend);
              glVertex3d(x[n1], y[n1], z[n1]);
              glColor4f(0, 0, 0+offb, blend);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
          }
    }
    glEnd();
  }
}
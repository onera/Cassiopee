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
# include "connector.h"

using namespace std;
using namespace K_FLD;

//==============================================================================
/* Met a jour les coordonnees des points frontieres d un maillage surfacique 
   Specifique a la creation de surfaces de projection en centres etendues
   issues de BCWall splittees*/
//==============================================================================
PyObject* K_CONNECTOR::modifyBorders(PyObject* self, PyObject* args)
{
  PyObject *a1;
  E_Int iminL0, jminL0, imaxL0, jmaxL0;
  if (!PYPARSETUPLE_(args, O_ IIII_,
                    &a1, &iminL0, &imaxL0, &jminL0, &jmaxL0))
  {
      return NULL;
  }

  // Check array of subzone of 1st layers of centers in extended centers mesh
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray(a1, varString1, f1, im1, jm1, km1, 
                                     cn1, eltType1, true); 
  if (res1 != 1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "modifyBorders: argument must define a structured array.");
    RELEASESHAREDB(res1, a1, f1, cn1); return NULL;
  }
  if (km1 != 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "modifyBorders: structured array must be 2D (k=1).");
    RELEASESHAREDB(res1, a1, f1, cn1); return NULL;
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1==-1 || posy1==-1 || posz1==-1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "modifyBorders: 1st arg must contain coordinates.");
    RELEASESHAREDB(res1, a1, f1, cn1); return NULL;
  }
  posx1++; posy1++; posz1++;
 
  E_Int npts = f1->getSize();
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", im1, jm1, km1);
  E_Float* fp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f(npts, 3, fp, true);
  f.setOneField(*f1,posx1,1);
  f.setOneField(*f1,posy1,2);
  f.setOneField(*f1,posz1,3);

  E_Float* xp = f.begin(posx1);
  E_Float* yp = f.begin(posy1);
  E_Float* zp = f.begin(posz1);
  E_Float* xp1 = f1->begin(posx1);
  E_Float* yp1 = f1->begin(posy1);
  E_Float* zp1 = f1->begin(posz1);
  E_Int ind, inci, incj;

  FldArrayI tag1(npts); tag1.setAllValuesAtNull();
  FldArrayI tag2(npts); tag2.setAllValuesAtNull();
  E_Int* tagp1 = tag1.begin(); E_Int* tagp2 = tag2.begin();
  E_Int iminL = E_Int(iminL0);
  E_Int jminL = E_Int(jminL0);
  E_Int imaxL = E_Int(imaxL0);
  E_Int jmaxL = E_Int(jmaxL0);
  
  if ( iminL > 1)
  {
    for (E_Int j=0; j < jm1; j++)
    {
      ind = (iminL-1)+j*im1;
      tagp1[ind]=-1; 
    }
  }
  if ( imaxL < im1 )
  {
    for (E_Int j=0; j < jm1; j++)
    {
      ind = (imaxL-1)+j*im1;
      tagp1[ind]= 1;
    }
  }
  
  if ( jminL > 1)
  {
    for (E_Int i=0; i < im1; i++)
    {
      ind = i+(jminL-1)*im1;
      tagp2[ind]=-2; 
    }
  }
  if ( jmaxL < jm1 )
  {
    for (E_Int i=0; i < im1; i++)
    {
      ind = i+(jmaxL-1)*im1;
      tagp2[ind]= 2;
    }
  }
  
  // Parcours des fenetres en i
  for (E_Int j = 0; j < jm1; j++)
    for (E_Int i = 0; i < im1; i++)
    {
      ind = i+j*im1;
      inci = 0; incj = 0;
      if (tagp1[ind]==-1) inci = 1;
      else if (tagp1[ind]==1) inci = -1;
      if (tagp2[ind]==-2) incj = im1;
      else if (tagp2[ind]==2) incj = -im1;
      if ( inci != 0 && incj == 0 )
      {
        xp[ind]=0.5*(xp1[ind]+xp1[ind+inci]);
        yp[ind]=0.5*(yp1[ind]+yp1[ind+inci]);
        zp[ind]=0.5*(zp1[ind]+zp1[ind+inci]);        
      }
      else if ( inci == 0 && incj != 0 )
      {
        xp[ind]=0.5*(xp1[ind]+xp1[ind+incj]);
        yp[ind]=0.5*(yp1[ind]+yp1[ind+incj]);
        zp[ind]=0.5*(zp1[ind]+zp1[ind+incj]);    
      }
      else if ( inci != 0 && incj != 0 )
      {
        xp[ind]=0.25*(xp1[ind]+xp1[ind+incj]+xp1[ind+inci]+xp1[ind+inci+incj]);
        yp[ind]=0.25*(yp1[ind]+yp1[ind+incj]+yp1[ind+inci]+yp1[ind+inci+incj]);
        zp[ind]=0.25*(zp1[ind]+zp1[ind+incj]+zp1[ind+inci]+zp1[ind+inci+incj]);
      }
    }
  
  RELEASESHAREDB(res1, a1, f1, cn1);
  return tpl;
}

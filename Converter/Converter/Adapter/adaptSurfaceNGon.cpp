/*    
    Copyright 2013-2023 Onera.

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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Adapt a surface NGon */
/* if INPUT is a NGON=bars, NFACE=polygon (Type A) -> OUTPUT is NGON=polygon, NFACE=NULL (Type B)
   if INPUT is NGON=polygon, NFACE=NULL (Type B) -> OUTPUT a NGON=bars, NFACE=polygon (Type A) */
//=============================================================================
PyObject* K_CONVERTER::adaptSurfaceNGon(PyObject* self, PyObject* args)
{
  PyObject* o; 
  if (!PYPARSETUPLE_(args, O_, &o)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(o, varString, f, ni, nj, nk, c, eltType);
  
  if (ret <= 0)
  { PyErr_SetString(PyExc_TypeError, "adaptSurfaceNGon: only for NGons."); return NULL; }

  if (ret == 1)
  { 
    PyErr_SetString(PyExc_TypeError, "adaptSurfaceNGon: only for NGons."); 
    RELEASESHAREDS(o, f);
    return NULL;
  }

  // Analyse input
  E_Int nelts = c->getNElts();
  PyObject* tpl = NULL;
  if (nelts == 0) // INPUT is type B  
  {
    E_Int nfacesB = c->getNFaces();

    E_Int neltsA = nfacesB;

  }
  else // INPUT is type A
  {
    E_Int isNGon = c->isNGon(); 
    E_Int nfacesB = nelts;
    E_Int neltsB = 0;
    E_Int sizeNGon = 0;
    E_Int sizeNFace = 0;
    tpl = K_ARRAY::buildArray3(f->getNfld(), varString, f->getSize(), neltsB, nfacesB, 
                               "NGON", sizeNGon, sizeNFace, isNGon, false, 3);
    K_FLD::FldArrayF* fo; K_FLD::FldArrayI* co;
    K_ARRAY::getFromArray3(tpl, fo, co);
    
    

    RELEASESHAREDU(tpl, fo, co);
  }

  RELEASESHAREDU(o, f, c);

  return tpl;
}

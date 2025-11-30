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

# include "kcore.h"
# include "Data.h"
# include "cplot.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Cree une zone structuree pour CPlot a partir d'un FldArray 
   IN: structF: coord + fields
   IN: varString: chaine des variables
   IN: posx, posy, posz: position des coords dans structF
   IN: ni,nj,nk: nbre de noeuds
   IN: zoneName: nom de la zone
   IN: zoneTags: tags de la zone (peut etre NULL)
   IN: referenceXX: sert de reference pour les variables
   IN: mustComplete: dit qu'il faut completer les variables (si appeler de add ou replace)
*/
//=============================================================================
StructZone* Data::createStructZone(FldArrayF* structF, char* varString,
                                   E_Int posx, E_Int posy, E_Int posz,
                                   E_Int ni, E_Int nj, E_Int nk,
                                   char* zoneName, char* zoneTags,
                                   E_Int referenceNfield, char** referenceVarNames,
                                   E_Int mustComplete)
{
  StructZone* sz = new StructZone(ptrState, createZoneImpl());
  StructZone& z = *sz;
  strcpy(z.zoneName, zoneName);

  z.ni = ni; z.nj = nj; z.nk = nk;
  z.npts = z.ni * z.nj * z.nk;
  if (referenceNfield != -1) z.nfield = referenceNfield;
  else z.nfield = structF->getNfld()-3;

  if (z.nk != 1) z.dim = 3;
  else if (z.nj != 1) z.dim = 2;
  else z.dim = 1;
  z.x = new E_Float[z.npts];
  memcpy(z.x, structF->begin(posx), z.npts*sizeof(E_Float));
  z.y = new E_Float[z.npts];
  memcpy(z.y, structF->begin(posy), z.npts*sizeof(E_Float));
  z.z = new E_Float[z.npts];
  memcpy(z.z, structF->begin(posz), z.npts*sizeof(E_Float));

  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int varsSize = vars.size();

  // Allocation of var fields
  reallocNFieldArrays(z.nfield);
  z.f = new double* [z.nfield];
  for (E_Int i = 0; i < z.nfield; i++) z.f[i] = NULL;
  z.varnames = new char* [z.nfield];
  for (E_Int n = 0; n < z.nfield; n++) 
    z.varnames[n] = new char [MAXSTRINGLENGTH];
  z.minf = new double [z.nfield];
  z.maxf = new double [z.nfield];

  if (referenceNfield != -1)
  {
    E_Int nall = 0;
    for (E_Int n = 0; n < referenceNfield; n++)
    {
      for (E_Int p = 0; p < varsSize; p++)
      {
        if (K_STRING::cmp(vars[p], referenceVarNames[n]) == 0)
        {
          nall++;
          z.f[n] = new E_Float[z.npts];
          memcpy(z.f[n], structF->begin(p+1), z.npts*sizeof(E_Float));
          strcpy(z.varnames[n], vars[p]); break;
        }
      }
    }
    nall = varsSize+referenceNfield-nall-3;
    for (E_Int n = 0; n < referenceNfield; n++)
    {
      if (z.f[n] == NULL)
      {
        z.f[n] = new E_Float[z.npts];
        if (K_STRING::cmp(referenceVarNames[n], "cellN") == 0)
        { for (int i = 0; i < z.npts; i++) z.f[n][i] = 1.; }
        else { for (int i = 0; i < z.npts; i++) z.f[n][i] = 0.; }
        strcpy(z.varnames[n], referenceVarNames[n]);
      }
    }
    // Complete all zones
    //printf("nall %d %d\n", nall, referenceNfield);
    if (nall > referenceNfield && mustComplete == 1)
    {
      // reallocate (previous zones)
      reallocNFieldArrays(nall);
      for (E_Int nz = 0; nz < _numberOfZones; nz++)
      {
        Zone* zp = _zones[nz];
        double** t = new double* [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t[i] = zp->f[i];
        delete [] zp->f; zp->f = t;

        char** t1 = new char* [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t1[i] = zp->varnames[i];
        for (E_Int i = zp->nfield; i < nall; i++) t1[i] = new char [MAXSTRINGLENGTH];
        delete [] zp->varnames; zp->varnames = t1;

        double* t2 = new double [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t2[i] = zp->minf[i];
        delete [] zp->minf; zp->minf = t2;
        
        double* t3 = new double [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t3[i] = zp->maxf[i];
        delete [] zp->maxf; zp->maxf = t3;
        
        zp->nfield = nall;
      }
      // reallocate current zone
      {
        double** t = new double* [nall];
        for (E_Int i = 0; i < z.nfield; i++) t[i] = z.f[i];
        delete [] z.f; z.f = t;

        char** t1 = new char* [nall];
        for (E_Int i = 0; i < z.nfield; i++) t1[i] = z.varnames[i];
        for (E_Int i = z.nfield; i < nall; i++) t1[i] = new char [MAXSTRINGLENGTH];
        delete [] z.varnames; z.varnames = t1;

        double* t2 = new double [nall];
        for (E_Int i = 0; i < z.nfield; i++) t2[i] = z.minf[i];
        delete [] z.minf; z.minf = t2;
        
        double* t3 = new double [nall];
        for (E_Int i = 0; i < z.nfield; i++) t3[i] = z.maxf[i];
        delete [] z.maxf; z.maxf = t3;

        z.nfield = nall;
      }
      
      nall = referenceNfield;
      for (E_Int p = 0; p < varsSize; p++)
      {
        E_Bool found = false; 
        for (E_Int n = 0; n < referenceNfield; n++)
        {
          if (K_STRING::cmp(vars[p], referenceVarNames[n]) == 0)
          {
            found = true; break;
          }
        }
        if (K_STRING::cmp(vars[p], "x") == 0 || K_STRING::cmp(vars[p], "y") == 0 ||
            K_STRING::cmp(vars[p], "z") == 0 || K_STRING::cmp(vars[p], "CoordinateX") == 0 ||
            K_STRING::cmp(vars[p], "CoordinateY") == 0 || K_STRING::cmp(vars[p], "CoordinateZ") == 0)
          found = true;
        
        if (found == false)
        {
          // ajoute vars[p] a la fin pour toutes les zones
          for (E_Int nz = 0; nz < _numberOfZones; nz++)
          {
            //printf("adding var %s (%d) in zone %d\n", vars[p], nall, nz);
            Zone* zp = _zones[nz];
            zp->f[nall] = new E_Float[zp->npts];
            if (K_STRING::cmp(vars[p], "cellN") == 0)
            { for (int i = 0; i < zp->npts; i++) zp->f[nall][i] = 1.; }
            else { for (int i = 0; i < zp->npts; i++) zp->f[nall][i] = 0.; }
            strcpy(zp->varnames[nall], vars[p]);
          }
          z.f[nall] = new E_Float[z.npts];
          memcpy(z.f[nall], structF->begin(p+1), z.npts*sizeof(E_Float));
          strcpy(z.varnames[nall], vars[p]);
          nall++;
        }
      }
    }
  }
  else // no reference zone
  {
    E_Int p = 0;
    for (E_Int n = 0; n < structF->getNfld(); n++)
    {
      if (n != posx-1 && n != posy-1 && n != posz-1)
      {
        z.f[p] = new E_Float[z.npts];
        memcpy(z.f[p], structF->begin(n+1), z.npts*sizeof(E_Float));
        strcpy(z.varnames[p], vars[n]); p++;
      }
    }
  }

  for (size_t n = 0; n < vars.size(); n++) delete [] vars[n];

  if (zoneTags != NULL)
  {
    strcpy(z.renderTag, zoneTags);
    codeFromRenderTag(z, z.renderTag, z.colorR, z.colorG, z.colorB, 
                      z.material, z.blending, z.meshOverlay, 
                      z.meshColorR, z.meshColorG, z.meshColorB, z.meshWidth,
                      z.shaderParam1, z.shaderParam2);
  }
  else
  {
    strcpy(z.renderTag, "None:None:None:None");
    z.colorR = -1; z.colorG = -1; z.colorB = -1; z.material = -1;
    z.blending = -1.; z.meshOverlay = 0;
    z.meshColorR = -1.; z.meshColorG = -1.; z.meshColorB = -1.; z.meshWidth = -1.; 
    z.shaderParam1 = 1.; z.shaderParam2 = 1.;
  }

  // Calcul les normales
  z.compNorm();

  z.activePlane = 0;
  z.iPlane = -1;
  z.jPlane = -1;
  z.kPlane = -1;
  z.iLine = 0;
  z.jLine = 0;
  z.kLine = 0;
  z.blank = 0;
  z.active = 1;
  z.selected = 0;
  findMinMax(&z); findFMinMax(&z);

  /* met les pointeurs sur _u_,_v_,_w_ pour le rendu texture */
  z.texu = NULL; z.texv = NULL; z.texw = NULL;
  z.regtexu = NULL; z.regtexv = NULL;
  if (z.material == 14)
  {
    for (E_Int n = 0; n < z.nfield; n++)
    { 
      if (strcmp(z.varnames[n], "texu") == 0 || strcmp(z.varnames[n], "_u_") == 0 || strcmp(z.varnames[n], "_U_") == 0) { z.texu = z.f[n]; }
      else if (strcmp(z.varnames[n], "texv") == 0 || strcmp(z.varnames[n], "_v_") == 0 || strcmp(z.varnames[n], "_V_") == 0) { z.texv = z.f[n]; }
      else if (strcmp(z.varnames[n], "texw") == 0 || strcmp(z.varnames[n], "_w_") == 0 || strcmp(z.varnames[n], "_W_") == 0) { z.texw = z.f[n]; }
    }
    if (z.texu == NULL)
    {
      // Create uniform texture field
      z.regtexu = new E_Float [z.npts];
      E_Float di = 1./(z.ni-1);
      for (E_Int j = 0; j < z.nj; j++)
        for (E_Int i = 0; i < z.ni; i++)
            z.regtexu[i+j*ni] = i*di;
      z.texu = z.regtexu;
    }
    if (z.texv == NULL)
    {
      // Create uniform texture field
      z.regtexv = new E_Float [z.npts];
      E_Float dj = 1./(z.nj-1);
      for (E_Int j = 0; j < z.nj; j++)
        for (E_Int i = 0; i < z.ni; i++)
            z.regtexv[i+j*ni] = j*dj;
      z.texv = z.regtexv;
    }
    
    if (z.texu != NULL && z.texv == NULL) z.texu = NULL;
    if (z.texu != NULL && z.texw == NULL) z.texw = z.texu;
  }

  return sz;
}

//=============================================================================
/* Cree une zone non-structuree pour CPlot a partir d'un FldArray 
   IN: unstrF: coord + fields
   IN: varString: chaine des variables
   IN: posx, posy, posz: position des coords dans unstrF
   IN: cn: connectivite
   IN: eltType: type d'element
   IN: zoneName: nom de la zone
   IN: zoneTags: tags de la zone (peut etre NULL)
*/
//=============================================================================
UnstructZone* Data::createUnstrZone(FldArrayF* unstrF, char* varString,
                                    E_Int posx, E_Int posy, E_Int posz,
                                    FldArrayI* cn, char* eltType,
                                    char* zoneName, char* zoneTags,
                                    E_Int referenceNfield, char** referenceVarNames,
                                    E_Int mustComplete)
{
  UnstructZone* uz = new UnstructZone(ptrState, createZoneImpl());
  UnstructZone& z = *uz;
  strcpy(z.zoneName, zoneName);
  z.npts = unstrF->getSize();
  z.np = z.npts;
  if (referenceNfield != -1) z.nfield = referenceNfield;
  else z.nfield = unstrF->getNfld()-3;

# if defined(__SHADERS__)
  const char* high_order_types[] = {
    "TRI_6", "TRI_9", "TRI_10", "TRI_12", "TRI_15", 
    "QUAD_8", "QUAD_9", "QUAD_12", "QUAD_16", "QUAD_P4_16", "QUAD_25",
    NULL
  };
  const short nb_nodes_per_elts[] = {
    6, 9, 10, 12, 15, 8, 9, 12, 16, 16, 25, -1  
  };
  unsigned short ind_type = 0;
  const char* pt_type = high_order_types[ind_type];
  while ( (pt_type != NULL) and (K_STRING::cmp(eltType, pt_type) != 0) )
  {
    ind_type++;
    pt_type = high_order_types[ind_type];
  }
  bool is_high_order = (pt_type != NULL);
  z._is_high_order = is_high_order;
# endif
  z.x = new E_Float[z.npts];
  memcpy(z.x, unstrF->begin(posx), z.npts*sizeof(E_Float));
  z.y = new E_Float[z.npts];
  memcpy(z.y, unstrF->begin(posy), z.npts*sizeof(E_Float));
  z.z = new E_Float[z.npts];
  memcpy(z.z, unstrF->begin(posz), z.npts*sizeof(E_Float));
  z.ne = cn->getSize();
  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int varsSize = vars.size();

  // Allocation of var fields
  reallocNFieldArrays(z.nfield);
  z.f = new double* [z.nfield];
  for (E_Int i = 0; i < z.nfield; i++) z.f[i] = NULL;
  z.varnames = new char* [z.nfield];
  for (E_Int n = 0; n < z.nfield; n++)
    z.varnames[n] = new char [MAXSTRINGLENGTH];
  z.minf = new double [z.nfield];
  z.maxf = new double [z.nfield];

  if (referenceNfield != -1)
  {
    E_Int nall = 0;
    for (E_Int n = 0; n < referenceNfield; n++)
    {
      for (E_Int p = 0; p < varsSize; p++)
      {
        if (K_STRING::cmp(vars[p], referenceVarNames[n]) == 0)
        {
          nall++;
          z.f[n] = new E_Float[z.npts];
          memcpy(z.f[n], unstrF->begin(p+1), z.npts*sizeof(E_Float));
          strcpy(z.varnames[n], vars[p]); break;
        }
      }
    }
    nall = varsSize+referenceNfield-nall-3;
    for (E_Int n = 0; n < referenceNfield; n++)
    {
      if (z.f[n] == NULL)
      {
        z.f[n] = new E_Float[z.npts];
        if (K_STRING::cmp(referenceVarNames[n], "cellN") == 0)
        { for (E_Int i = 0; i < z.npts; i++) z.f[n][i] = 1.; }
        else { for (E_Int i = 0; i < z.npts; i++) z.f[n][i] = 0.; }
        strcpy(z.varnames[n], referenceVarNames[n]);
      }
    }
    // Complete all zones
    if (nall > referenceNfield && mustComplete == 1)
    {
      // reallocate (previous zones)
      reallocNFieldArrays(nall);
      for (E_Int nz = 0; nz < _numberOfZones; nz++)
      {
        Zone* zp = _zones[nz];
        double** t = new double* [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t[i] = zp->f[i];
        delete [] zp->f; zp->f = t;

        char** t1 = new char* [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t1[i] = zp->varnames[i];
        for (E_Int i = zp->nfield; i < nall; i++) t1[i] = new char [MAXSTRINGLENGTH];
        delete [] zp->varnames; zp->varnames = t1;

        double* t2 = new double [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t2[i] = zp->minf[i];
        delete [] zp->minf; zp->minf = t2;
        
        double* t3 = new double [nall];
        for (E_Int i = 0; i < zp->nfield; i++) t3[i] = zp->maxf[i];
        delete [] zp->maxf; zp->maxf = t3;
        
        zp->nfield = nall;
      }
      // reallocate current zone
      {
        double** t = new double* [nall];
        for (E_Int i = 0; i < z.nfield; i++) t[i] = z.f[i];
        delete [] z.f; z.f = t;

        char** t1 = new char* [nall];
        for (E_Int i = 0; i < z.nfield; i++) t1[i] = z.varnames[i];
        for (E_Int i = z.nfield; i < nall; i++) t1[i] = new char [MAXSTRINGLENGTH];
        delete [] z.varnames; z.varnames = t1;

        double* t2 = new double [nall];
        for (E_Int i = 0; i < z.nfield; i++) t2[i] = z.minf[i];
        delete [] z.minf; z.minf = t2;
        
        double* t3 = new double [nall];
        for (E_Int i = 0; i < z.nfield; i++) t3[i] = z.maxf[i];
        delete [] z.maxf; z.maxf = t3;

        z.nfield = nall;
      }
      
      nall = referenceNfield;
      for (E_Int p = 0; p < varsSize; p++)
      {
        E_Bool found = false; 
        for (E_Int n = 0; n < referenceNfield; n++)
        {
          if (K_STRING::cmp(vars[p], referenceVarNames[n]) == 0)
          {
            found = true; break;
          }
        }
        if (K_STRING::cmp(vars[p], "x") == 0 || K_STRING::cmp(vars[p], "y") == 0 ||
            K_STRING::cmp(vars[p], "z") == 0 || K_STRING::cmp(vars[p], "CoordinateX") == 0 ||
            K_STRING::cmp(vars[p], "CoordinateY") == 0 || K_STRING::cmp(vars[p], "CoordinateZ") == 0)
          found = true;
            
        if (found == false)
        {
          // ajoute vars[p] a la fin pour toutes les zones
          for (E_Int nz = 0; nz < _numberOfZones; nz++)
          {
            //printf("adding var %s (%d) in zone %d\n", vars[p], nall, nz);
            Zone* zp = _zones[nz];
            zp->f[nall] = new E_Float[zp->npts];
            if (K_STRING::cmp(vars[p], "cellN") == 0)
            { for (E_Int i = 0; i < zp->npts; i++) zp->f[nall][i] = 1.; }
            else { for (E_Int i = 0; i < zp->npts; i++) zp->f[nall][i] = 0.; }
            strcpy(zp->varnames[nall], vars[p]);
          }
          z.f[nall] = new E_Float[z.npts];
          memcpy(z.f[nall], unstrF->begin(p+1), z.npts*sizeof(E_Float));
          strcpy(z.varnames[nall], vars[p]);
          nall++;
        }
      }
    }
  }
  else // no reference zone
  {
    E_Int p = 0;
    for (E_Int n = 0; n < unstrF->getNfld(); n++)
    {
      if (n != posx-1 && n != posy-1 && n != posz-1)      
      {
        z.f[p] = new E_Float[z.npts];
        memcpy(z.f[p], unstrF->begin(n+1), z.npts*sizeof(E_Float));
        strcpy(z.varnames[p], vars[n]); p++;
      }
    }
  }

  for (E_Int n = 0; n < varsSize; n++) delete [] vars[n];
  
  if (zoneTags != NULL)
  {
    strcpy(z.renderTag, zoneTags);
    codeFromRenderTag(z, z.renderTag, z.colorR, z.colorG, z.colorB, 
                      z.material, z.blending, z.meshOverlay, 
                      z.meshColorR, z.meshColorG, z.meshColorB, z.meshWidth,
                      z.shaderParam1, z.shaderParam2);
  }
  else
  {
    strcpy(z.renderTag, "None:None:None:None:None:None:None:None");
    z.colorR = -1.; z.colorG = -1.; z.colorB = -1.; z.material = -1;
    z.blending = -1.; z.meshOverlay = 0; 
    z.meshColorR = -1.; z.meshColorG = -1.; z.meshColorB = -1.; z.meshWidth = -1.;
    z.shaderParam1 = 1.; z.shaderParam2 = 1.;
  }

  /* Explore connectivities */
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (size_t i = 0; i < eltTypes.size(); i++)
  {
    if (K_STRING::cmp(eltTypes[i], "NODE") == 0)
    {
      z.eltType.push_back(0);
      z.eltSize.push_back(0);
      z.dim = 0;
    }
    else if (K_STRING::cmp(eltTypes[i], "BAR") == 0)
    {
      z.eltType.push_back(1);
      z.eltSize.push_back(2);
      z.dim = 1;
    }
    else if (K_STRING::cmp(eltTypes[i], "TRI") == 0)
    {
      z.eltType.push_back(2);
      z.eltSize.push_back(3);
      z.dim = 2;
    }
    else if (K_STRING::cmp(eltTypes[i], "QUAD") == 0)
    {
      z.eltType.push_back(3);
      z.eltSize.push_back(4);
      z.dim = 2;
    }
    else if (K_STRING::cmp(eltTypes[i], "TETRA") == 0)
    {
      z.eltType.push_back(4);
      z.eltSize.push_back(4);
      z.dim = 3;
    }
    else if (K_STRING::cmp(eltTypes[i], "PENTA") == 0)
    {
      z.eltType.push_back(5);
      z.eltSize.push_back(6);
      z.dim = 3;
    }
    else if (K_STRING::cmp(eltTypes[i], "PYRA") == 0)
    {
      z.eltType.push_back(6);
      z.eltSize.push_back(5);
      z.dim = 3;
    }
    else if (K_STRING::cmp(eltTypes[i], "HEXA") == 0)
    {
      z.eltType.push_back(7);
      z.eltSize.push_back(8);
      z.dim = 3;
    }
    else if (K_STRING::cmp(eltTypes[i], "NGON") == 0)
    {
      z.ne = cn->getNElts();
      z.nec.push_back(z.ne);
      z.eltType.push_back(10);
      z.eltSize.push_back(1);
      z.dim = 3;
    }
    else if (not z._is_high_order)
    {
      printf("Warning: element type is unknown. Set to TRI.\n");
      z.eltType.push_back(2);
      z.eltSize.push_back(3);
      z.dim = 2;
    }
  }
# if defined(__SHADERS__)
  if (is_high_order)
  {
    z.eltSize.push_back(nb_nodes_per_elts[ind_type]);
    z.eltType.push_back((ind_type < 5 ? 2 : 3));
    z.dim = 2;
  }
# endif

  for (size_t i = 0; i < eltTypes.size(); i++) delete [] eltTypes[i];
  
  // copie par access universel
  if (z.eltType[0] == 10) // NGON
  {
    E_Int nfaces = cn->getNFaces();
    E_Int nelts = cn->getNElts();
    E_Int size1 = cn->getSizeNGon();
    E_Int size2 = cn->getSizeNFace();
    if (cn->getNGonType() == 3)
    {
      size1 += nfaces; size2 += nelts;
    }
    z.connect.push_back(new E_Int[size1+size2+4]);
    E_Int* zconnect = z.connect[0];
    zconnect[0] = nfaces;
    zconnect[1] = size1;
    zconnect[size1+2] = nelts;
    zconnect[size1+3] = size2;
    E_Int* znp = zconnect+2;
    E_Int* fp; E_Int size;
    E_Int* ngon = cn->getNGon();
    E_Int* nface = cn->getNFace();
    E_Int* indPG = cn->getIndPG();
    E_Int* indPH = cn->getIndPH();
    for (E_Int i = 0; i < nfaces; i++)
    {
      fp = cn->getFace(i, size, ngon, indPG);
      znp[0] = size;
      for (E_Int v = 0; v < size; v++) znp[v+1] = fp[v];
      znp += size+1;
    }
    znp += 2;
    for (E_Int i = 0; i < nelts; i++)
    {
      fp = cn->getElt(i, size, nface, indPH);
      znp[0] = size;
      for (E_Int v = 0; v < size; v++) znp[v+1] = std::abs(fp[v]);
      znp += size+1;
    }
  }
  else // BE
  {
    E_Int ncon = cn->getNConnect();
    E_Int neTot = 0;
    for (E_Int nc = 0; nc < ncon; nc++)
    {
      FldArrayI& cm = *(cn->getConnect(nc));
      E_Int nvpe = cm.getNfld();
      E_Int nelts = cm.getSize();
      neTot += nelts;
      E_Int size = nelts * nvpe;
      z.connect.push_back(new E_Int[size]);
      E_Int* znp = z.connect[nc];
      for (E_Int n = 0; n < nvpe; n++)
        for (E_Int i = 0; i < nelts; i++)
          znp[i+n*nelts] = cm(i, n+1);
      z.nec.push_back(nelts);
    }
    z.ne = neTot;
  }


  z.posFaces = NULL;

  if (z.eltType[0] == 10) // NGONS
  {
    // calcul posFaces (position des faces dans connect)
    E_Int* zconnect = z.connect[0];
    E_Int nfaces = NFACES(zconnect);
    z.posFaces = new E_Int[nfaces];
    E_Int c = POSFACES(zconnect); E_Int l;
    for (E_Int i = 0; i < nfaces; i++)
    {
      z.posFaces[i] = c; l = zconnect[c]; c += l+1;
    }

    // calcul le nombre d'elements 1D et 2D
    E_Int nelts = NELTS(zconnect);
    
    c = POSELTS(zconnect);
    E_Int dim, s, c1, c2;
    z.nelts1D = 0; z.nelts2D = 0;
    for (E_Int i = 0; i < nelts; i++)
    {
      l = zconnect[c]; // nbre de faces
      dim = 0;
      for (E_Int j = 0; j < l; j++)
      {
        s = z.posFaces[zconnect[c+j+1]-1];
        dim = max(dim, zconnect[s]);
      }
      if (dim == 1) { z.nelts1D++; }
      else if (dim == 2) { z.nelts2D++; }
      c += l+1;
    }
    
    //printf("1D: %d, 2D: %d\n", z.nelts1D, z.nelts2D);
    if (z.nelts1D > 0) z.posElts1D = new E_Int[z.nelts1D];
    if (z.nelts2D > 0) z.posElts2D = new E_Int[z.nelts2D];
    c = POSELTS(zconnect);
    c1 = 0; c2 = 0;
    for (E_Int i = 0; i < nelts; i++)
    {
      l = zconnect[c]; // nbre de faces
      dim = 0;
      for (E_Int j = 0; j < l; j++)
      {
        s = z.posFaces[zconnect[c+j+1]-1];
        dim = max(dim, zconnect[s]);
      }
      if (dim == 1) { z.posElts1D[c1] = c; c1++; }
      else if (dim == 2) { z.posElts2D[c2] = c; c2++; }
      c += l+1;
    }
  }

  z.compNorm();
  z.blank = 0;
  z.active = 1;
  z.selected = 0;
  findMinMax(&z); findFMinMax(&z);

  /* met les pointeurs sur u,v,w pour le rendu texture */
  z.texu = NULL; z.texv = NULL; z.texw = NULL;
  z.regtexu = NULL; z.regtexv = NULL;
  if (z.material == 14)
  {
    for (E_Int n = 0; n < z.nfield; n++)
    { 
      if (strcmp(z.varnames[n], "texu") == 0 || strcmp(z.varnames[n], "_u_") == 0 || strcmp(z.varnames[n], "_U_") == 0) { z.texu = z.f[n]; }
      else if (strcmp(z.varnames[n], "texv") == 0 || strcmp(z.varnames[n], "_v_") == 0 || strcmp(z.varnames[n], "_V_") == 0) { z.texv = z.f[n]; }
      else if (strcmp(z.varnames[n], "texw") == 0 || strcmp(z.varnames[n], "_w_") == 0 || strcmp(z.varnames[n], "_W_") == 0) { z.texw = z.f[n]; }
    }
    
    if (z.texu != NULL && z.texv == NULL) z.texu = NULL;
    if (z.texu != NULL && z.texw == NULL) z.texw = z.texu;
  }
  return uz;
}

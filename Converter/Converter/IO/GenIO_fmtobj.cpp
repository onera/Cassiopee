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

// Formated Obj file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* objread */
// Manque: commentaires
// Nettoyage de f a la fin (vertex reellement utilises)
//=============================================================================
E_Int K_IO::GenIO::objread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t;
  E_Int ti, ti1, ti2, i, j, nv, nf, nt, nq, k;
  char buf[256];
  FldArrayF* f;
  FldArrayI* cn_q=NULL;
  FldArrayI* cn_t=NULL;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: objread: cannot open file %s.\n", file);
    return 1;
  }

  // Recherche du nombre de materiaux differents
  vector<char*> matnames;
  E_Int nm = 0; E_Int nom = 0; res = 1;
  vector<E_Int> matindir; matindir.reserve(200);
  char buffer[256];
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "USEMTL ");
    if (res == 1) 
    { 
      readWord(ptrFile, buffer);
      // exist?
      bool exist = false;
      for (size_t i = 0; i < matnames.size(); i++)
      {
        if (strcmp(buffer, matnames[i]) == 0) { exist=true; matindir.push_back(i); break; }
      }
      if (exist == false)
      {
        char* name = new char [256];
        strcpy(name, buffer);
        matnames.push_back(name);
        matindir.push_back(nm);
        nm++;
      }
      nom++;
    }
  }
  //printf("I found %d materials (occurences=%d)\n", nm, nom);
  //for (E_Int i = 0; i < nm; i++) printf("material %d: %s\n", i, matnames[i]);
  //for (E_Int i = 0; i < nom; i++) printf("occurence %d: %I64d\n", i, matpos[i]); 

  // Recherche du nombre de V
  KFSEEK(ptrFile, 0, SEEK_SET);
  nv = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "V "); nv++;
  }
  nv = nv-1;
  //printf("I found %d vertices\n", nv);

  // Recherche si normales
  KFSEEK(ptrFile, 0, SEEK_SET);
  res = readGivenKeyword(ptrFile, "VN ");
  //bool normalPresent = false;
  //if (res == 1) normalPresent = true;
  //if (normalPresent) printf("I found normals\n");
  //else printf("I found no normal\n");

  // Recherche du nombre de VT
  KFSEEK(ptrFile, 0, SEEK_SET);
  E_Int nvt = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "VT "); nvt++;
  }
  nvt = nvt-1;
  bool uvPresent = false;
  if (nvt > 0) { uvPresent = true; }
  //printf("I found %d texcoord VT\n", nvt);

  // DBX -> force no uv
  //uvPresent = false;

  // Lecture des faces
  KFSEEK(ptrFile, 0, SEEK_SET);
  nf = 0; res = 1;
  while (res == 1)
  {
    res = readGivenKeyword(ptrFile, "F "); nf++;
  }
  nf = nf-1;
  //printf("I found %d faces\n", nf);

  // Type des faces
  KFSEEK(ptrFile, 0, SEEK_SET);
  FldArrayI elts(nf);
  res = 1;
  res = readGivenKeyword(ptrFile, "F ");
  i = 0; nt = 0; nq = 0;
  while (res >= 1)
  {
    for (j = 0; j < 4; j++)
    {
      res = readWord(ptrFile, buf); 
    }
    if (res >= 1 && strcmp(buf, "f") == 0)
    {elts[i] = 3; nt++; }
    else if ((res >= 1 || res == 0) && (buf[0] < 48 || buf[0] > 57))
    {elts[i] = 3; nt++; res = readGivenKeyword(ptrFile, "F "); }
    else if (res == 0 && buf[0] >= 48 && buf[0] <= 57)
    {elts[i] = 4; nq++; res = readGivenKeyword(ptrFile, "F ");}
    else if (res == -1)
    {elts[i] = 3; nt++; res = readGivenKeyword(ptrFile, "F "); }
    else { elts[i] = 4; nq++; res = readGivenKeyword(ptrFile, "F ");}
    i++;
  }
  
  //printf("I found %d tri %d quads\n", nt, nq);
  
  FldArrayI fv_q; FldArrayI fv_t;
  FldArrayI fvt_q; FldArrayI fvt_t;
  FldArrayI fom_q; FldArrayI fom_t; // occurence materiau

  // Look for quads
  if (nq > 0)
  {
    cn_q = new FldArrayI(nq, 4);
    fv_q.malloc(nq, 4); // vertices
    fvt_q.malloc(nq, 4); // vt
    fom_q.malloc(nq); // occurence materiau

    KFSEEK(ptrFile, 0, SEEK_SET);
    E_Int nom = -1;
    res = 1; i = 0; k = 0;
    res = readGivenKeyword(ptrFile, "F ", "USEMTL ");
    
    while (res >= 1)
    {
      if (res == 2) { nom++;}
      else if (res == 1)
      {
        if (elts[k] == 4)
        {
          fom_q[i] = nom;
          for (j = 1; j <= 4; j++)
          {
            res = readIntTuple3(ptrFile, ti, ti1, ti2);
            fv_q(i,j) = ti; fvt_q(i,j) = ti1;
          }
          i++;
        }
        k++;
      }
      res = readGivenKeyword(ptrFile, "F ", "USEMTL ");
    }
  }

  // Look for tri
  if (nt > 0)
  {
    cn_t = new FldArrayI(nt, 3);
    fv_t.malloc(nt, 3);
    fvt_t.malloc(nt, 3);
    fom_t.malloc(nt);

    KFSEEK(ptrFile, 0, SEEK_SET);
    E_Int nom = -1;
    res = 1; i = 0; k = 0;
    res = readGivenKeyword(ptrFile, "F ", "USEMTL ");
    while (res >= 1)
    {
      if (res == 2) { nom++;}
      else if (res == 1)
      {
        if (elts[k] == 3)
        {
          fom_t[i] = nom;
          for (j = 1; j <= 3; j++)
          {
            res = readIntTuple3(ptrFile, ti, ti1, ti2);
            fv_t(i,j) = ti; fvt_t(i,j) = ti1;
          }
          //printf("%d %d %d\n", (*cn_t)(i,1), (*cn_t)(i,2), (*cn_t)(i,3));
          i++;
        }
        k++;
      }
      res = readGivenKeyword(ptrFile, "F ", "USEMTL ");
    }
  }

  // Traitement specifique pour les indices negatifs (ne marche pas si par bloc)
  for (E_Int i = 0; i < nq; i++)
  {
    for (E_Int j = 1; j <= 4; j++)
    {
      if (fv_q(i,j) < 1) fv_q(i,j) = nv+fv_q(i,j);
      if (fvt_q(i,j) < 1) fvt_q(i,j) = nvt+fvt_q(i,j);
    }
  }
  for (E_Int i = 0; i < nt; i++)
  {
    for (E_Int j = 1; j <= 3; j++)
    {
      if (fv_t(i,j) < 1) fv_t(i,j) = nv+fv_t(i,j);
      if (fvt_t(i,j) < 1) fvt_t(i,j) = nvt+fvt_t(i,j);
    }
  }

  // map v -> vt
  FldArrayI map(nv*20);
  map.setAllValuesAtNull();
  FldArrayI nmap(nv);
  nmap.setAllValuesAtNull();

  E_Int v, vt;
  bool exist;

  //printf("Computing maps\n");
  for (E_Int i = 0; i < nq; i++)
  {
    for (E_Int j = 1; j <= 4; j++)
    {
      v = fv_q(i,j)-1; vt = fvt_q(i,j);
      
      if (vt > 0)
      {
        // existe deja dans la map?
        exist = false;
        for (E_Int k = 0; k < nmap[v]; k++)
        {
          if (map[v + k*nv] == vt) { exist=true; break; }
        }
        if (exist == false && nmap[v] < 19) { map[v + nmap[v]*nv] = vt; nmap[v]++; } 
      }
    }
  }

  for (E_Int i = 0; i < nt; i++)
  {
    for (E_Int j = 1; j <= 3; j++)
    {
      v = fv_t(i,j)-1; vt = fvt_t(i,j);
      if (vt > 0)
      {
        // existe deja dans la map?
        exist = false;
        for (E_Int k = 0; k < nmap[v]; k++)
        {
          if (map[v + k*nv] == vt) { exist=true; break; }
        }
        if (exist == false && nmap[v] < 19) { map[v + nmap[v]*nv] = vt; nmap[v]++; } 
      }
    }
  }

  // indir v -> vnew
  //printf("Computing indir\n");
  FldArrayF indir(nv);
  indir[0] = 0;
  for (E_Int i = 1; i < nv; i++)
    indir[i] = indir[i-1]+max(nmap[i-1],1);

  // mise a plat et dimensionnement de f
  E_Int size = 0;
  for (E_Int i = 0; i < nv; i++) size += max(nmap[i],1);
  //printf("full size=%d\n", size);

  if (uvPresent) f = new FldArrayF(size, 5);
  else f = new FldArrayF(nv, 3);

  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);

  // Lecture vertices
  FldArrayF* coord = NULL;
  E_Float* px = NULL;
  E_Float* py = NULL;
  E_Float* pz = NULL;
  if (uvPresent)
  {
    coord = new FldArrayF(nv, 3);
    px = coord->begin(1);
    py = coord->begin(2);
    pz = coord->begin(3);
  }
  else
  {
    px = fx; py = fy; pz = fz;
  }

  KFSEEK(ptrFile, 0, SEEK_SET);
  i = 0; res = 1;
  res = readGivenKeyword(ptrFile, "V ");
  while (res == 1)
  {
    res = readDouble(ptrFile, t, -1); px[i] = t; //printf("%f ", t);
    res = readDouble(ptrFile, t, -1); py[i] = t; //printf("%f ", t);
    res = readDouble(ptrFile, t, -1); pz[i] = t; i++; //printf("%f\n", t);
    //if (res == 0) res = 1; else res = 0;
    res = readGivenKeyword(ptrFile, "V ");
  }

  // Lecture VT -> used as Vertices u,v
  FldArrayF texcoord;
  if (uvPresent)
  {
    texcoord.malloc(nvt,2);
    E_Float* tu = texcoord.begin(1);
    E_Float* tv = texcoord.begin(2);

    KFSEEK(ptrFile, 0, SEEK_SET);
    i = 0;
    res = readGivenKeyword(ptrFile, "VT ");
    while (res == 1)
    {
      res = readDouble(ptrFile, t, -1); tu[i] = t; //printf("%d %d: %f ", i, nv, t);
      res = readDouble(ptrFile, t, -1); tv[i] = t; i++; //printf("%f\n", t);
      res = readGivenKeyword(ptrFile, "VT ");
    }
    // Traitement specifique vt (tiles)
    E_Int ip;
    for (E_Int i = 0; i < nvt; i++)
    {
      if (tu[i] > 1) { ip = E_Int(tu[i]); tu[i] = tu[i]-ip; }
      else if (tu[i] < 0) { ip = E_Int(tu[i]); tu[i] = tu[i]-ip+1; }
      if (tv[i] > 1) { ip = E_Int(tv[i]); tv[i] = tv[i]-ip; }
      else if (tv[i] < 0) { ip = E_Int(tv[i]); tv[i] = tv[i]-ip+1; }
    }  
  }

  // mise a plat de f
  if (uvPresent)
  {
    E_Float* tu = f->begin(4);
    E_Float* tv = f->begin(5);
    E_Int ind;
    E_Int pt = 0;
    E_Float* pu = texcoord.begin(1);
    E_Float* pv = texcoord.begin(2); 
    for (E_Int i = 0; i < nv; i++)
    {
      if (nmap[i] == 0) { fx[pt] = px[i]; fy[pt] = py[i]; fz[pt] = pz[i]; tu[pt] = 0.; tv[pt] = 0.; pt++; }
      for (E_Int j = 0; j < nmap[i]; j++)
      {
        ind = map[i+j*nv]-1;
        fx[pt] = px[i]; fy[pt] = py[i]; fz[pt] = pz[i]; tu[pt] = pu[ind]; tv[pt] = pv[ind]; pt++;
      }
    } 
  }

  // nouvelles connectivites
  for (E_Int i = 0; i < nq; i++)
  {
    for (E_Int j = 1; j <= 4; j++)
    {
      v = fv_q(i,j)-1;
      vt = fvt_q(i,j);
      E_Int k = 0;
      for (k = 0; k < nmap[v]; k++)
      {
        if (map[v + k*nv] == vt) break;
      }
      if (k == nmap[v]) k = 0; // not found
      (*cn_q)(i,j) = indir[v]+k+1;
    }
  }

  for (E_Int i = 0; i < nt; i++)
  {
    for (E_Int j = 1; j <= 3; j++)
    {
      v = fv_t(i,j)-1;
      vt = fvt_t(i,j);
      E_Int k = 0;
      for (k = 0; k < nmap[v]; k++)
      {
        if (map[v + k*nv] == vt) break;
      }
      if (k == nmap[v]) k = 0; // not found
      (*cn_t)(i,j) = indir[v]+k+1;
    }
  }

  nmap.malloc(0);
  map.malloc(0);

  // sortie une zone par materiau (TRI)
  if (nt > 0)
  {
    // attribue chaque face a un materiau
    FldArrayI matface(nt);
    matface.setAllValuesAt(-1);

    for (E_Int i = 0; i < nt; i++)
    {
      if (fom_t[i] == -1) matface[i] = -1;
      else matface[i] = matindir[fom_t[i]];
    }

    for (E_Int m = -1; m < nm; m++)
    {
      // Compte les faces
      E_Int nf = 0;
      for (E_Int i = 0; i < nt; i++)
      {
        if (matface[i] == m) nf++; 
      }
      //printf("mat=%d, TRI nf=%d\n", m, nf);
      if (nf > 0)
      {
        // Dimensionne
        FldArrayI* cn = new FldArrayI(nf, 3);
        // Rempli
        E_Int ff = 0;
        for (E_Int i = 0; i < nt; i++)
        {
          if (matface[i] == m) 
          {
            for (E_Int j = 1; j <= 3; j++) (*cn)(ff,j) = (*cn_t)(i,j);
            //printf("%d %d %d\n", (*cn)(f,1), (*cn)(ff,2), (*cn)(ff,3));
            ff++;
          }
        }
      
        // clean unreferenced vertices in f
        FldArrayF* f2 = new FldArrayF();
        FldArrayI* cn2 = new FldArrayI();
        K_CONNECT::cleanUnreferencedVertexBasic(*f, *cn, *f2, *cn2);
        delete cn;

        unstructField.push_back(f2);
        eltType.push_back(2);
        connect.push_back(cn2);
        char* zoneName = new char [128];
        if (m == -1) strcpy(zoneName, "NoMatTRI");
        else sprintf(zoneName, "%sTRI", matnames[m]);
        zoneNames.push_back(zoneName);
      }
    }
    delete cn_t;
  }

  if (nq > 0)
  {
    // attribue chaque face a un materiau
    FldArrayI matface(nq);
    matface.setAllValuesAt(-1);

    for (E_Int i = 0; i < nq; i++)
    {
      if (fom_q[i] == -1) matface[i] = -1;
      else matface[i] = matindir[fom_q[i]];
    }
    
    for (E_Int m = -1; m < nm; m++)
    {
      // Compte les faces
      E_Int nf = 0;
      for (E_Int i = 0; i < nq; i++)
      {
        if (matface[i] == m) nf++; 
      }
      //printf("mat=%d, QUAD nf=%d\n", m, nf);
      if (nf > 0)
      {
        // Dimensionne
        FldArrayI* cn = new FldArrayI(nf, 4);
        // Rempli
        E_Int ff = 0;
        for (E_Int i = 0; i < nq; i++)
        {
          if (matface[i] == m) 
          {
            for (E_Int j = 1; j <= 4; j++) (*cn)(ff,j) = (*cn_q)(i,j);
            ff++;
          } 
        }
        // push back
        // clean unreferenced vertices
        FldArrayF* f2 = new FldArrayF();
        FldArrayI* cn2 = new FldArrayI();
        K_CONNECT::cleanUnreferencedVertexBasic(*f, *cn, *f2, *cn2);
        delete cn;

        unstructField.push_back(f2);
        eltType.push_back(3);
        connect.push_back(cn2);
        char* zoneName = new char [128];
        if (m == -1) strcpy(zoneName, "NoMatQUAD");
        else sprintf(zoneName, "%sQUAD", matnames[m]);
        zoneNames.push_back(zoneName);
      }
    }
    delete cn_q;
  }

  delete f;
  for (size_t i = 0; i < matnames.size(); i++) delete [] matnames[i];

  varString = new char [16];
  if (uvPresent == false)
  { strcpy(varString, "x,y,z"); }
  else { strcpy(varString, "x,y,z,u,v"); }
  fclose(ptrFile);
  
  return 0;
}

//=============================================================================
// Only write triangles, quads, tetra, hexa meshes. Others are discarded.
//=============================================================================
E_Int K_IO::GenIO::objwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  E_Int zone;
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (zone = 0; zone < nzone; zone++)
  {
    // triangles, quads supported
    if (eltType[zone] == 2 || eltType[zone] == 3)
      nvalidZones++;
    else
      printf("Warning: objwrite: zone %d not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: objwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Ouverture fichier
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: objwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Build writing data format
  char format1[124]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt);
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format of data
  sprintf(format1,"v %s%s%s\n", dataFmt, dataFmt, dataFmtl);

  // Write zones
  E_Int ng = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    // vertices
    FldArrayF& field = *unstructField[i];
    E_Float* fx = field.begin(posx);
    E_Float* fy = field.begin(posy);
    E_Float* fz = field.begin(posz);
    for (E_Int n = 0; n < field.getSize(); n++)
    {
      fprintf(ptrFile, format1, fx[n], fy[n], fz[n]);
    }
    if (eltType[i] == 2)
    {
      // faces
      FldArrayI& cn = *connect[i];
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        fprintf(ptrFile, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", 
                cn1[n]+ng, cn1[n]+ng, cn1[n]+ng,
                cn2[n]+ng, cn2[n]+ng, cn2[n]+ng,
                cn3[n]+ng, cn3[n]+ng, cn3[n]+ng);
      }
    }
    if (eltType[i] == 3)
    {
      // faces
      FldArrayI& cn = *connect[i];
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      E_Int* cn4 = cn.begin(4);
      for (E_Int n = 0; n < cn.getSize(); n++)
      {
        fprintf(ptrFile, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n", 
                cn1[n]+ng, cn1[n]+ng, cn1[n]+ng,
                cn2[n]+ng, cn2[n]+ng, cn2[n]+ng,
                cn3[n]+ng, cn3[n]+ng, cn3[n]+ng,
                cn4[n]+ng, cn4[n]+ng, cn4[n]+ng);
      }
    }
    ng += field.getSize();
  }

  fclose(ptrFile);
  return 0;
}

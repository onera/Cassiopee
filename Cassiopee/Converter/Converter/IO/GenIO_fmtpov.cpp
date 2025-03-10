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

// Formated Pov (Povray) file support

# include <stdio.h>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6normunstructsurf_(E_Int& nbt, E_Int& sizecoord, E_Int* cn, 
                           E_Float* coordx, E_Float* coordy, E_Float* coordz,
                           E_Float* nsurf);
}

//=============================================================================
/* povread */
//=============================================================================
E_Int K_IO::GenIO::povread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t;
  E_Int size;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: povread: cannot open file %s.\n", file);
    return 1;
  }

  // Boucles sur les domaines
  res = 1; 

  while (res != 0)
  {
    // Recherche le mot cle mesh2 
    res = readGivenKeyword(ptrFile, "MESH2");
    if (res == 0) goto end;
    res = readGivenKeyword(ptrFile, "{");
    if (res == 0) goto end;

    // Coordonnees
    res = readGivenKeyword(ptrFile, "VERTEX_VECTORS");
    if (res == 0) goto end;
    res = readGivenKeyword(ptrFile, "{");
    if (res == 0) goto end;
    res = readDouble(ptrFile, t, -1);
    size = E_Int(t);
    //printf("size = " SF_D_ "\n", size);
    FldArrayF* field = new FldArrayF(size,3);
    FldArrayF& f = *field;
    
    for (E_Int i = 0; i < size; i++)
    {
      res = readDouble(ptrFile, t, -1); f(i,1) = t;
      res = readDouble(ptrFile, t, -1); f(i,2) = t;
      res = readDouble(ptrFile, t, -1); f(i,3) = t;
      //printf(SF_F3_ "\n", f(i,1), f(i,2), f(i,3));
    }
    res = readGivenKeyword(ptrFile, "}");
    if (res == 0) {delete field; goto end;}
    
    // Connectivite
    res = readGivenKeyword(ptrFile, "FACE_INDICES");
    if (res == 0) {delete field; goto end;}
    res = readGivenKeyword(ptrFile, "{");
    if (res == 0) {delete field; goto end;}
    res = readDouble(ptrFile, t, -1);
    size = E_Int(t);
    //printf("size = " SF_D_ "\n", size);
    FldArrayI* cn = new FldArrayI(size, 3);
    FldArrayI& c = *cn;
    
    for (E_Int i = 0; i < size; i++)
    {
      res = readDouble(ptrFile, t, -1); c(i,1) = int(t+1);
      res = readDouble(ptrFile, t, -1); c(i,2) = int(t+1);
      res = readDouble(ptrFile, t, -1); c(i,3) = int(t+1);
      //printf(" SF_D3_ "\n", c(i,1), c(i,2), c(i,3));
    }

    unstructField.push_back(field);
    connect.push_back(cn);
    eltType.push_back(2);

    res = readGivenKeyword(ptrFile, "}");
    if (res == 0) {delete cn; delete field; goto end;}
    res = readGivenKeyword(ptrFile, "}");
    if (res == 0) {delete cn; delete field; goto end;}
  }
  
  end:
  // Cree les noms des zones
  for (size_t i = 0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }
  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Only write triangles arrays. Others are discarded.
// Write it in mesh2 format.
//=============================================================================
E_Int K_IO::GenIO::povwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes, 
  vector<char*>& zoneNames, E_Int colormap)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltTypes[zone][0] == 2) // triangles 
      nvalidZones++;
    else
      printf("Warning: povwrite: zone " SF_D_ " not written (not a triangle zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  E_Int posd = 1;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: povwrite: zones do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  if (colormap > 0)
  {
    posd = colormap;
    if (posd > unstructField[0]->getNfld())
    {
      printf("Warning: povwrite: zones do not have this field. Only coordinates are written.");
      colormap = 0;
    }
  }

  char format1[126], format2[128];
  char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format
  sprintf(format1,"<%s,%s,%s>,", dataFmt, dataFmt, dataFmtl);
  sprintf(format2,"<%s,%s,%s>\n}\n", dataFmt, dataFmt, dataFmtl);
  // Find min max of first variable if possible
  E_Float fmin = K_CONST::E_MAX_FLOAT;
  E_Float fmax = -K_CONST::E_MAX_FLOAT;
  E_Float delta = 0.;
  if (colormap > 0)
  {
    for (E_Int zone = 0; zone < nzone; zone++)
    {
      FldArrayF& f = *unstructField[zone];
      E_Float* fp = f.begin(posd);
      for (E_Int l = 0; l < f.getSize(); l++)
      {
        fmin = K_FUNC::E_min(fmin, fp[l]);
        fmax = K_FUNC::E_max(fmax, fp[l]);
      }
    }
  }
  //printf(SF_F2_ \n", fmin, fmax);

  // Create colormap if necessary
  FldArrayF rgb;
  E_Int col1, col2, col3;
  E_Float N = 1.;
  if (colormap > 0)
  {
    createColormap(colormap, rgb);
    N = rgb.getSize();
    E_Float dd = K_FUNC::E_max(fmax-fmin, 1.e-13);
    delta = N/dd;
    //printf("N = " SF_F2_ " \n", N, delta);
  }

  // Open file
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: povwrite: I can't open file %s.\n", file);
    return 1;
  }

  // Mesh name
  char meshName[256];
  { 
    E_Int m = 0;
    while (file[m] != '.' && file[m] != '\0')
    {meshName[m] = file[m]; m++;}
    meshName[m] = '\0';
  }

  for (E_Int zone = 0; zone < nzone; zone++)
  {
    E_Int nv = unstructField[zone]->getSize();
    FldArrayF& f = *unstructField[zone];
    E_Float* fx = f.begin(posx);
    E_Float* fy = f.begin(posy);
    E_Float* fz = f.begin(posz);

    // Vertices
    fprintf(ptrFile, "#declare %s" SF_D_ "=mesh2 {\n", meshName, zone);
    fprintf(ptrFile, "    vertex_vectors {\n");
    fprintf(ptrFile, "       " SF_D_ ",\n", nv);
    E_Int l = 0;
    for (E_Int i = 0; i < nv-1; i++)
    {
      fprintf(ptrFile, format1, fx[i], fy[i], fz[i]);
      if (l > 3) 
      {
        l = 0; fprintf(ptrFile, "\n");
      }
      l++;
    }
    fprintf(ptrFile, format2, fx[nv-1], fy[nv-1], fz[nv-1]);

    // Normals
    E_Float nx, ny, nz, nt;
    E_Int ne = connect[zone]->getSize();
    FldArrayI& c = *connect[zone];
    FldArrayF normals(ne, 3);
    k6normunstructsurf_(ne, nv, c.begin(), 
                        f.begin(posx), f.begin(posy), f.begin(posz),
                        normals.begin());
    E_Float* npx = normals.begin(1);
    E_Float* npy = normals.begin(2);
    E_Float* npz = normals.begin(3);
    for (E_Int i = 0; i <  ne; i++)
    {
      nx = npx[i]; ny = npy[i]; nz = npz[i];
      nt = 1./K_FUNC::E_max(sqrt(nx*nx + ny*ny + nz*nz), 1.e-12);
      npx[i] = nx * nt;
      npy[i] = ny * nt;
      npz[i] = nz * nt;
    }
    
    vector< vector<E_Int> > cVE(nv);
    K_CONNECT::connectEV2VE(c, cVE);
    FldArrayF n(nv, 3);
    E_Float* n1 = n.begin(1);
    E_Float* n2 = n.begin(2);
    E_Float* n3 = n.begin(3);
    for (E_Int i = 0; i < nv; i++)
    {
      vector<E_Int>& cVEi = cVE[i];
      E_Int sizecVE = K_FUNC::E_max(cVEi.size(), E_Int(1));
      E_Int cVESize = cVEi.size();
      n1[i] = 0; n2[i] = 0; n3[i] = 0;
      for (E_Int j = 0; j < cVESize; j++)
      {
        n1[i] = n1[i] + npx[cVEi[j]];
        n2[i] = n2[i] + npy[cVEi[j]];
        n3[i] = n3[i] + npz[cVEi[j]];
      }
      n1[i] = n1[i] / sizecVE;
      n2[i] = n2[i] / sizecVE;
      n3[i] = n3[i] / sizecVE;
    }
    fprintf(ptrFile, "    normal_vectors {\n");
    fprintf(ptrFile, "       " SF_D_ ",\n", nv);
    
    l = 0;
    for (E_Int i = 0; i < nv-1; i++)
    {
      fprintf(ptrFile, format1, n1[i], n2[i], n3[i]);
      if (l > 3) 
      {
        l = 0; 
        fprintf(ptrFile, "\n");
      }
      l++;
    }
    fprintf(ptrFile, format2, 
            n1[nv-1], n2[nv-1], n3[nv-1]);

    // Texture list
    if (colormap > 0)
    {
      fprintf(ptrFile, "    texture_list {\n");
      fprintf(ptrFile, "       " SF_D_ ",\n", rgb.getSize());
      for (E_Int i = 0; i < rgb.getSize(); i++)
        fprintf(ptrFile, "     texture{pigment{rgb<" SF_F_ "," SF_F_ "," SF_F_ ">}}\n",
                rgb(i,1), rgb(i,2), rgb(i,3));
      fprintf(ptrFile, " }\n");
    }

    // Connectivity
    fprintf(ptrFile, "    face_indices {\n");
    fprintf(ptrFile, "       " SF_D_ ",\n", ne);
    l = 0;
    if (colormap == 0)
    {
      for (E_Int i = 0; i < ne-1; i++)
      {
        fprintf(ptrFile, "<" SF_D_ "," SF_D_ "," SF_D_ ">,",
                c(i,1)-1, c(i,2)-1, c(i,3)-1);
        if (l > 3) 
        {
          l = 0; 
          fprintf(ptrFile, "\n");
        }
        l++;
      }
      fprintf(ptrFile, "<" SF_D_ "," SF_D_ "," SF_D_ ">\n}\n", 
              c(ne-1,1)-1, c(ne-1,2)-1, c(ne-1,3)-1);
    }
    else
    {
      for (E_Int i = 0; i < ne-1; i++)
      {
        col1 = E_Int(K_FUNC::E_min((f(c(i,1)-1, posd) - fmin)*delta, N-1)); 
        col2 = E_Int(K_FUNC::E_min((f(c(i,2)-1, posd) - fmin)*delta, N-1)); 
        col3 = E_Int(K_FUNC::E_min((f(c(i,3)-1, posd) - fmin)*delta, N-1));
        fprintf(ptrFile, "<" SF_D_ "," SF_D_ "," SF_D_ ">," SF_D_ "," SF_D_ "," SF_D_ ",", 
                c(i,1)-1, c(i,2)-1, c(i,3)-1,
                col1, col2, col3);
        if (l > 3) 
        {
          l = 0; 
          fprintf(ptrFile, "\n");
        }
        l++;
      }
      col1 = E_Int(K_FUNC::E_min((f(c(ne-1,1)-1, posd) - fmin)*delta, N-1)); 
      col2 = E_Int(K_FUNC::E_min((f(c(ne-1,2)-1, posd) - fmin)*delta, N-1)); 
      col3 = E_Int(K_FUNC::E_min((f(c(ne-1,3)-1, posd) - fmin)*delta, N-1)); 
      fprintf(ptrFile, "<" SF_D_ "," SF_D_ "," SF_D_ ">," SF_D_ "," SF_D_ "," SF_D_ "\n}\n", 
              c(ne-1,1)-1, c(ne-1,2)-1, c(ne-1,3)-1,
              col1, col2, col3);
    }

    fprintf(ptrFile, "    }\n"); // end of this zone
  }  

  if (nzone > 1)
  {
    fprintf(ptrFile, "#declare %s=union {\n", meshName);
    for (E_Int zone = 0; zone < nzone; zone++)
    {
      fprintf(ptrFile, "object {%s" SF_D_ "}\n", meshName, zone);
    }
    fprintf(ptrFile, "}\n");
  }
  else
    fprintf(ptrFile, "#declare %s=object {%s0}\n", meshName, meshName);

  fclose(ptrFile);
  return 0;
}

//=============================================================================
void K_IO::GenIO::createColormap(E_Int colormap, FldArrayF& rgb)
{
  // blue to red colormap
  {
    E_Int N = 20;
    rgb.malloc(2*N,3);
    E_Int c = 0;
    for (E_Int r = 0; r < N; r++)
    {
      rgb(c,1) = 0.;
      rgb(c,2) = r*1./(N-1);
      rgb(c,3) = 1. - r*1./(N-1);
      c++;
    }
    for (E_Int r = 0; r < N; r++)
    {
      rgb(c,1) = r*1./(N-1);
      rgb(c,2) = 1.-r*1./(N-1);
      rgb(c,3) = 0.; 
      c++;
    }
  }
}

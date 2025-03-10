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

// Binary PLY (Stanford) file support

# include "GenIO.h"
# include "Array/Array.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   plyread 
*/
//=============================================================================
E_Int K_IO::GenIO::plyread( 
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: plyread: cannot open file %s.\n", file);
    return 1;
  }
  
  // Header
  char dummy[81];
  fscanf(ptrFile, "%s", dummy);
  if (strcmp(dummy, "ply") != 0) { fclose(ptrFile); return 1; }
  
  fscanf(ptrFile, "%s", dummy);
  if (strcmp(dummy, "format") != 0) { fclose(ptrFile); return 1; }

  fscanf(ptrFile, "%s", dummy);
  if (strcmp(dummy, "ascii") == 0) { fclose(ptrFile); return 1; }
  int endian = 0;
  if (strcmp(dummy, "binary_little_endian") == 0) endian = 0;
  if (strcmp(dummy, "binary_big_endian") == 0) endian = 1;
  E_Int machineEndian = machineEndianess();
  E_Int convertEndian = 0;
  if (machineEndian != endian) convertEndian = 1;
  //if (convertEndian == 1) printf("Need to translate endians\n");

  // Recherche element
  E_Int found = 0;
  while (found == 0)
  { 
    fscanf(ptrFile, "%s", dummy);
    if (strcmp(dummy, "element") == 0) found = 1; 
  }
  fscanf(ptrFile, "%s", dummy); // type element
  
  if (strcmp(dummy, "vertex") != 0) { fclose(ptrFile); return 1; }
  fscanf(ptrFile, "%s", dummy); // nbre de vertex
  E_Int nd = atoi(dummy);
  
  // Recherche property
  E_Int typeVar[20];
  char nameVar[20][128];
  E_Int nvar, type;

  fscanf(ptrFile, "%s", dummy);
  nvar = 0;
  while (strcmp(dummy, "property") == 0 && nvar < 20)
  {
    fscanf(ptrFile, "%s", dummy); // type
    type = 0;
    if (strcmp(dummy, "float") == 0) type = 0;
    else if (strcmp(dummy, "double") == 0) type = 1;
    typeVar[nvar] = type;
    fscanf(ptrFile, "%s", dummy); // variable
    strcpy(nameVar[nvar], dummy);
    nvar++;
    fscanf(ptrFile, "%s", dummy);
  }
  
  // recherche element face
  if (strcmp(dummy, "element") == 0) found = 1;
  else found = 0;

  while (found == 0)
  { 
    fscanf(ptrFile, "%s", dummy);
    if (strcmp(dummy, "element") == 0) found = 1; 
  }
  fscanf(ptrFile, "%s", dummy); // type element
  if (strcmp(dummy, "face") != 0) {fclose(ptrFile); return 1;}

  fscanf(ptrFile, "%s", dummy); // nbre de face
  E_Int ne = atoi(dummy);
  
  // Recherche property
  //E_Int typec = 0;
  fscanf(ptrFile, "%s", dummy);
  if (strcmp(dummy, "property") != 0) {fclose(ptrFile); return 1;}
  fscanf(ptrFile, "%s", dummy);
  if (strcmp(dummy, "list") != 0) {fclose(ptrFile); return 1;}
  fscanf(ptrFile, "%s", dummy); // type
  //if (strcmp(dummy, "uchar") == 0) typec = 0;
  //else if (strcmp(dummy, "int") == 0) typec = 1;
  //else if (strcmp(dummy, "uint") == 0) typec = 2;
  fscanf(ptrFile, "%s", dummy); // variable
  
  while(strcmp(dummy, "end_header") != 0)
  {
    fscanf(ptrFile, "%s", dummy);
  }
  char buf[1];
  fread(buf, sizeof(char), 1, ptrFile); // rid of \n

  // Champ
  FldArrayF* f = new FldArrayF(nd, 3);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);
  
  //for (E_Int i = 0; i < nvar; i++) printf("%s\n", nameVar[i]);

  float fbuf[1]; double dbuf[1];
  if (convertEndian == 0)
  {
    for (E_Int i = 0; i < nd; i++)
    {
      for (E_Int nv = 0; nv < nvar; nv++)
      {
        if (typeVar[nv] == 0 && strcmp(nameVar[nv], "x") == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); fx[i] = fbuf[0];
        }
        else if (typeVar[nv] == 0 && strcmp(nameVar[nv], "y") == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); fy[i] = fbuf[0];
        }
        else if (typeVar[nv] == 0 && strcmp(nameVar[nv], "z") == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); fz[i] = fbuf[0];
        }
        else if (typeVar[nv] == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); // lost
        }
        else if (typeVar[nv] == 1 && strcmp(nameVar[nv], "x") == 0)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); fx[i] = dbuf[0];
        }
        else if (typeVar[nv] == 1 && strcmp(nameVar[nv], "y") == 0)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); fy[i] = dbuf[0];
        }
        else if (typeVar[nv] == 1 && strcmp(nameVar[nv], "z") == 0)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); fz[i] = dbuf[0];
        }
        else if (typeVar[nv] == 1)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); // lost
        }
        //else printf("missed\n");
      }
    }
  }
  else // endian conversion
  {
    for (E_Int i = 0; i < nd; i++)
    {
      for (E_Int nv = 0; nv < nvar; nv++)
      {
        if (typeVar[nv] == 0 && strcmp(nameVar[nv], "x") == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); fx[i] = FBE(fbuf[0]);
        }
        else if (typeVar[nv] == 0 && strcmp(nameVar[nv], "y") == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); fy[i] = FBE(fbuf[0]);
        }
        else if (typeVar[nv] == 0 && strcmp(nameVar[nv], "z") == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); fz[i] = FBE(fbuf[0]);
        }
        else if (typeVar[nv] == 0)
        {
          fread(fbuf, sizeof(float), 1, ptrFile); // lost
        }
        else if (typeVar[nv] == 1 && strcmp(nameVar[nv], "x") == 0)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); fx[i] = DBE(dbuf[0]);
        }
        else if (typeVar[nv] == 1 && strcmp(nameVar[nv], "y") == 0)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); fy[i] = DBE(dbuf[0]);
        }
        else if (typeVar[nv] == 1 && strcmp(nameVar[nv], "z") == 0)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); fz[i] = DBE(dbuf[0]);
        }
        else if (typeVar[nv] == 1)
        {
          fread(dbuf, sizeof(double), 1, ptrFile); // lost
        }
        //else printf("missed\n");
      }
    }
  }
  
  // Lecture Connectivite
  // Connectivite TRI interlace avec les QUADS
  FldArrayI* ct = new FldArrayI(ne, 3);
  FldArrayI* cq = new FldArrayI(ne, 4);

  E_Int* ct1 = ct->begin(1);
  E_Int* ct2 = ct->begin(2);
  E_Int* ct3 = ct->begin(3);
  E_Int* cq1 = cq->begin(1);
  E_Int* cq2 = cq->begin(2);
  E_Int* cq3 = cq->begin(3);
  E_Int* cq4 = cq->begin(4);

  // Je ne lis que uchar uint pour l'instant
  unsigned char buf1[1];
  unsigned int buf2[5];
  int nt = 0;
  int ntri = 0; int nquad = 0;

  for (E_Int i = 0; i < ne; i++)
  {
    fread(buf1, sizeof(unsigned char), 1, ptrFile);
    nt = buf1[0];
    if (nt == 3) // TRI
    {
      fread(buf2, sizeof(unsigned int), 3, ptrFile);
      ct1[ntri] = buf2[0]+1;
      ct2[ntri] = buf2[1]+1;
      ct3[ntri] = buf2[2]+1; ntri++;
    }
    else if (nt == 4) // QUAD
    {
      fread(buf2, sizeof(unsigned int), 4, ptrFile);
      cq1[nquad] = buf2[0]+1;
      cq2[nquad] = buf2[1]+1;
      cq3[nquad] = buf2[2]+1;
      cq4[nquad] = buf2[3]+1; nquad++;
    }
    else fread(buf2, sizeof(unsigned int), nt, ptrFile); // lost 
  }

  if (ntri > 0 && nquad == 0)
  {
    unstructField.push_back(f);
    ct->reAllocMat(ntri, 3);
    connect.push_back(ct);
    eltType.push_back(2);
    K_CONNECT::cleanConnectivity(1, 2, 3, 
                                 1.e-14,  "TRI",
                                 *unstructField[0], *connect[0]);
  }
  else if (nquad > 0 && ntri == 0)
  {
    unstructField.push_back(f);
    cq->reAllocMat(nquad, 4);
    connect.push_back(cq);
    eltType.push_back(3);
    K_CONNECT::cleanConnectivity(1, 2, 3, 
                                 1.e-14,  "QUAD",
                                 *unstructField[0], *connect[0]);
  }
  else if (ntri == 0 && nquad == 0)
  {
    printf("Warning: plyread: file contains no element.\n");
  }
  else
  {
    unstructField.push_back(f);
    FldArrayF* f2 = new FldArrayF(f->getSize(), f->getNfld());
    *f2 = *f;
    ct->reAllocMat(ntri, 3);
    connect.push_back(ct);
    eltType.push_back(2);
    K_CONNECT::cleanConnectivity(1, 2, 3, 
                                 1.e-14,  "TRI",
                                 *unstructField[0], *connect[0]);
    
    unstructField.push_back(f2);
    cq->reAllocMat(nquad, 4);
    connect.push_back(cq);
    eltType.push_back(3);
    K_CONNECT::cleanConnectivity(1, 2, 3, 
                                 1.e-14,  "QUAD",
                                 *unstructField[1], *connect[1]);
  }

  // Cree les noms de zones
  E_Int unstructFieldSize = unstructField.size();
  for (E_Int i = 0; i < unstructFieldSize; i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, i);
    zoneNames.push_back(zoneName);
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  return 0;
}

//=============================================================================
// Ecrit 1 array TRI et 1 array QUAD (mixed)
//=============================================================================
E_Int K_IO::GenIO::plywrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int ntri = -1; E_Int nquad = -1;

  E_Int missed = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltTypes[zone][0] == 2) // triangles
    { if (ntri == -1) ntri = zone; else missed = 1; }
    if (eltTypes[zone][0] == 3) // Quads
    { if (nquad == -1) nquad = zone; else missed = 1; }
  }

  if (ntri == -1 && nquad == -1) return 1;
  if (missed == 1) printf("Warning: plywrite: only first TRI and QUAD zones are written.\n");
  
  // Zone must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: plywrite: zone do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  // Open file
  FILE* ptrFile = fopen(file, "wb");
  if (ptrFile == NULL) 
  {
    printf("Warning: plywrite: I can't open file %s.\n", file);
    return 1;
  }

  // Check machine endianess
  E_Int endianess = machineEndianess();

  // Header
  fprintf(ptrFile, "ply\n");
  // toujours little endian
  fprintf(ptrFile, "format binary_little_endian 1.0\n");
  fprintf(ptrFile, "comment Written by Cassiopee\n");

  // Tableaux
  FldArrayF* at = NULL; FldArrayI* ct = NULL;
  FldArrayF* aq = NULL; FldArrayI* cq = NULL;
  if (ntri != -1) { at = unstructField[ntri]; ct = connect[ntri]; }
  if (nquad != -1) { aq = unstructField[nquad]; cq = connect[nquad]; }
  
  // Vertex
  E_Int nv = 0;
  if (at != NULL) nv += at->getSize();
  if (aq != NULL) nv += aq->getSize();
  fprintf(ptrFile, "element vertex " SF_D_ "\n", nv);
  fprintf(ptrFile, "property double x\n");
  fprintf(ptrFile, "property double y\n");
  fprintf(ptrFile, "property double z\n");

  // Faces
  E_Int nf = 0;
  if (ct != NULL) nf += ct->getSize();
  if (cq != NULL) nf += cq->getSize();
  fprintf(ptrFile, "element face " SF_D_ "\n", nf);
  fprintf(ptrFile, "property list uchar uint vertex_index\n");
  
  fprintf(ptrFile, "end_header\n");

  // Ecriture vertex
  double dbuf[3];
  E_Int ndec = 0;
  if (at != NULL)
  {
    E_Float* fx = at->begin(posx);
    E_Float* fy = at->begin(posy);
    E_Float* fz = at->begin(posz);
    E_Int nd = at->getSize();
    ndec = nd;

    if (endianess == 0)
    {
      for (E_Int i = 0; i < nd; i++)
      {
        dbuf[0] = fx[i]; dbuf[1] = fy[i]; dbuf[2] = fz[i];
        fwrite(dbuf, sizeof(double), 3, ptrFile);
      }
    }
    else
    {
      for (E_Int i = 0; i < nd; i++)
      {
        dbuf[0] = DBE(fx[i]); dbuf[1] = DBE(fy[i]); dbuf[2] = DBE(fz[i]);
        fwrite(dbuf, sizeof(double), 3, ptrFile);
      }
    }
  }
  if (aq != NULL)
  {
    E_Float* fx = aq->begin(posx);
    E_Float* fy = aq->begin(posy);
    E_Float* fz = aq->begin(posz);
    E_Int nd = aq->getSize();

    if (endianess == 0)
    {
      for (E_Int i = 0; i < nd; i++)
      {
        dbuf[0] = fx[i]; dbuf[1] = fy[i]; dbuf[2] = fz[i];
        fwrite(dbuf, sizeof(double), 3, ptrFile);
      }
    }
    else
    {
      for (E_Int i = 0; i < nd; i++)
      {
        dbuf[0] = DBE(fx[i]); dbuf[1] = DBE(fy[i]); dbuf[2] = DBE(fz[i]);
        fwrite(dbuf, sizeof(double), 3, ptrFile);
      }
    }
  }

  // Ecriture connectivite
  unsigned char ubuf[1]; unsigned int ibuf[4];
  if (ct != NULL)
  {
    E_Int ne = ct->getSize();
    FldArrayI& c = *ct;
    if (endianess == 0)
    {
      ubuf[0] = 3;
      for (E_Int i = 0; i < ne; i++)
      {
        fwrite(ubuf, sizeof(unsigned char), 1, ptrFile);
        ibuf[0] = c(i,1)-1; ibuf[1] = c(i,2)-1; ibuf[2] = c(i,3)-1;
        fwrite(ibuf, sizeof(unsigned int), 3, ptrFile);
      }
    }
    else
    {
      E_Int nt = 3;
      ubuf[0] = IBE(nt); // not sure
      for (E_Int i = 0; i < ne; i++)
      {
        fwrite(ubuf, sizeof(unsigned char), 1, ptrFile);
        ibuf[0] = IBE(c(i,1)-1); ibuf[1] = IBE(c(i,2)-1); ibuf[2] = IBE(c(i,3)-1);
        fwrite(ibuf, sizeof(unsigned int), 3, ptrFile);
      }
    }
  }
  if (cq != NULL)
  {
    E_Int ne = cq->getSize();
    FldArrayI& c = *cq;
    if (endianess == 0)
    {
      ubuf[0] = 4;
      for (E_Int i = 0; i < ne; i++)
      {
        fwrite(ubuf, sizeof(unsigned char), 1, ptrFile);
        ibuf[0] = c(i,1)-1+ndec; ibuf[1] = c(i,2)-1+ndec; 
        ibuf[2] = c(i,3)-1+ndec; ibuf[3] = c(i,4)-1+ndec;
        fwrite(ibuf, sizeof(unsigned int), 4, ptrFile);
      }
    }
    else
    {
      E_Int nt = 4;
      ubuf[0] = IBE(nt); // not sure
      for (E_Int i = 0; i < ne; i++)
      {
        fwrite(ubuf, sizeof(unsigned char), 1, ptrFile);
        ibuf[0] = IBE(c(i,1)-1+ndec); ibuf[1] = IBE(c(i,2)-1+ndec); 
        ibuf[2] = IBE(c(i,3)-1+ndec); ibuf[3] = IBE(c(i,4)-1+ndec);
        fwrite(ibuf, sizeof(unsigned int), 4, ptrFile);
      }
    }
  }

  fclose(ptrFile);
  return 0;
}

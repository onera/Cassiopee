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

// Binary vtk (legacy) file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

// retourne 0 si ok, 1 sinon
int readLine(FILE* ptrFile, char* buf)
{
  char c = '1'; E_Int i = 0;
  while (c != '\n')
  { 
    if (feof(ptrFile)) return 1;
    c = fgetc(ptrFile);
    buf[i] = c; i++; 
  }
  if (i > 0) buf[i-1] = '\0';
  else buf[0] = '\0';
  //printf("LINE: %s\n", buf);
  return 0;
}

void getFirstWord(char* buf, char* word, E_Int& length)
{
  E_Int i = 0;
  while (buf[i] != '\n')
  {
    word[i] = buf[i]; i++;
  }
  if (i > 0) { word[i-1] = '\0'; length = i-1; }
  else { word[0] = '\0'; length = 0; }
  printf("WORD:%s\n", word);
}

//========================================================================
// read POINTS
void readPoints(FILE* ptrFile, E_Bool changeEndian, E_Bool formated,
                E_Int& npts, float*& bf, double*& bd)
{
  bf = NULL; bd = NULL;
  char buf[256]; char type[256];
  readLine(ptrFile, buf);
  E_Int l = strlen(buf);
  E_Int count = 0;
  while (buf[count] != ' ' && buf[count] != '\n') count++;
  for (E_Int i = 0; i <= l-count; i++) buf[i] = buf[i+count+1];
  
  count = 0;
  while (buf[count] != ' ' && buf[count] != '\n') count++;
  for (E_Int i = 0; i <= l-count; i++) type[i] = buf[i+count+1];
  buf[count] = '\0';

  npts = 0;
  sscanf(buf, SF_D_, &npts);
  printf("POINTS npts: " SF_D_ ", type=%s\n", npts, type);

  if (strcmp(type, "double") == 0)
  {
    bd = new double [npts*3];
    if (formated == false)
    {
        fread(bd, sizeof(double), 3*npts, ptrFile);
        fgetc(ptrFile); // skip \n
    }
    else
    {
        for (E_Int i = 0; i < npts; i++) 
            fscanf(ptrFile, "%lf %lf %lf\n", &bd[3*i], &bd[3*i+1], &bd[3*i+2]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < 3*npts; i++) bd[i] = DBE(bd[i]);
    }
    printf("point0 %lf %f %lf\n", bd[0], bd[1], bd[2]);
    printf("point1 %lf %lf %lf\n", bd[3], bd[4], bd[5]);
  }
  else
  {
    bf = new float [npts*3];
    if (formated == false)
    {
        fread(bf, sizeof(float), 3*npts, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < npts; i++) 
            fscanf(ptrFile, "%f %f %f\n", &bf[3*i], &bf[3*i+1], &bf[3*i+2]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < 3*npts; i++) bf[i] = FBE(bf[i]);
    }
    printf("point0 %f %f %f\n", bf[0], bf[1], bf[2]);
    printf("point1 %f %f %f\n", bf[3], bf[4], bf[5]);
  }
   
}

//===========================================================================
void readCells(FILE* ptrFile, E_Bool changeEndian, E_Bool formated, E_Int& ncells, E_Int& size, int*& cells)
{
  char buf[256]; char buf2[256]; char cellType[256];
  readLine(ptrFile, buf);
  
  // recupere le cell type
  E_Int i = 0;
  while (buf[i] != '\0' && buf[i] != ' ') { cellType[i] = buf[i]; i++; }
  cellType[i] = '\0';
  E_Int prev = strlen(cellType)+1;
  printf("cell Type=%s\n", cellType);

  // Recupere le nombre de cellules
  i = prev;
  while (buf[i] != '\0' && buf[i] != ' ') { buf2[i-prev] = buf[i]; i++; }
  buf2[i-prev] = '\0';
  ncells = 0;
  sscanf(buf2, SF_D_, &ncells);
  printf("ncells: " SF_D_ "\n", ncells);
  // recupere la taille totale de la connectivite
  i += 1;
  E_Int is = i;
  while (buf[i] != '\0' && buf[i] != ' ') { buf2[i-is] = buf[i]; i++; }
  buf2[i-is] = '\0';
  sscanf(buf2, SF_D_, &size);
  printf("size: " SF_D_ "\n", size);
  cells = new int [size];
  if (formated == false)
  {
    fread(cells, sizeof(int), size, ptrFile);
    fgetc(ptrFile); // avoid \n
  }
  else
  {
    // force to triangle for now
    //E_Int npic;
    //for (E_Int i = 0; i < size; i++) 
    //{
    //    fscanf(ptrFile, "%d %d %d %d\n", &npic, &cells);    
    //}
  }
  if (changeEndian)
  {
    for (E_Int i = 0; i < size; i++) cells[i] = IBE(cells[i]);
  }
  printf("cells=%d\n", cells[0]);
}

//===========================================================================
void readCellTypes(FILE* ptrFile, E_Bool changeEndian, E_Bool formated, 
                   E_Int ncells, int*& cellTypes)
{
  char buf[256];
  readLine(ptrFile, buf);
  cellTypes = new int [ncells];
  if (formated == false)
  {
    fread(cellTypes, sizeof(int), ncells, ptrFile);
    fgetc(ptrFile); // avoid \n
  }
  else
  {
    for (E_Int i = 0; i < ncells; i++) fscanf(ptrFile, "%d", &cellTypes[i]);
  }
  if (changeEndian)
  {
    for (E_Int i = 0; i < ncells; i++) cellTypes[i] = IBE(cellTypes[i]);
  } 
  printf("cellTypes=%d\n", cellTypes[0]);
}

//=========================================================================
void readType(FILE* ptrFile, char* rtype, E_Int& n1, E_Int& n2)
{
  char buf[256]; char buf2[256];
  readLine(ptrFile, buf);
  
  // recupere le mot cle
  E_Int i = 0;
  while (buf[i] != '\0' && buf[i] != ' ') { rtype[i] = buf[i]; i++; }
  rtype[i] = '\0';
  E_Int prev = strlen(rtype)+1;
  printf("rtype=%s\n", rtype);

  // Recupere le nombre de cellules
  i = prev;
  while (buf[i] != '\0' && buf[i] != ' ') { buf2[i-prev] = buf[i]; i++; }
  buf2[i-prev] = '\0';
  n1 = 0;
  sscanf(buf2, SF_D_, &n1);
  printf("n1: " SF_D_ "\n", n1);

  // recupere la taille totale de la connectivite
  i += 1;
  E_Int is = i;
  while (buf[i] != '\0' && buf[i] != ' ') { buf2[i-is] = buf[i]; i++; }
  buf2[i-is] = '\0';
  sscanf(buf2, SF_D_, &n2);
  printf("n2: " SF_D_ "\n", n2);
}

//===========================================================================
void readPolygons(FILE* ptrFile, E_Bool changeEndian, E_Bool formated,
                  E_Int& ncells, E_Int& size, int*& cells)
{
  char celltype[256];
  readType(ptrFile, celltype, ncells, size);
  printf("ncells = " SF_D_ ", size=" SF_D_ "\n", ncells, size);
  cells = new int [size];
  E_Int num;

  if (formated)
  {
      E_Int c = 0;
      for (E_Int i = 0; i < ncells; i++) 
      {
        fscanf(ptrFile, SF_D_ " ", &num);
        cells[c] = num; c++;
        if (num == 2) fscanf(ptrFile, "%d %d\n", &cells[c], &cells[c+1]);
        else if (num == 3) fscanf(ptrFile, "%d %d %d\n", &cells[c], &cells[c+1], &cells[c+2]);
        else if (num == 4) fscanf(ptrFile, "%d %d %d %d\n", &cells[c], &cells[c+1], &cells[c+2], &cells[c+3]);
        c += num;
      }
  }
  else
  {
      fread(cells, sizeof(E_Int), size, ptrFile);
  }
  printf("cells[0]=%d %d %d\n", cells[0], cells[1], cells[2]);
}

//========================================================================
void readScalar(FILE* ptrFile, E_Bool changeEndian, E_Bool formated,
                E_Int npts, char* varName, float*& fieldf, double*& fieldd)
{ 
  fieldf = NULL; fieldd = NULL;

  char buf[256];

  // recupere dataName, dataType, numComp
  char dataType[256]; E_Int numComp;
  fscanf(ptrFile, " %s %s " SF_D_ "\n", varName, dataType, &numComp);
  printf("varName=%s, dataType=%s, numCom=" SF_D_ "\n", varName, dataType, numComp);

  // Passe lookup table
  readLine(ptrFile, buf);

  if (strcmp(dataType, "float") == 0)
  {
    fieldf = new float [npts];
    if (formated == false)
    {
        fread(fieldf, sizeof(float), npts, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%f", &fieldf[i]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < npts; i++) fieldf[i] = FBE(fieldf[i]);
    }
  }
  else if (strcmp(dataType, "double") == 0)
  {
    fieldd = new double [npts];
    if (formated == false)
    {
        fread(fieldd, sizeof(double), npts, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%lf", &fieldd[i]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < npts; i++) fieldd[i] = DBE(fieldd[i]);
    }
  }
}

//==========================================================================
void readVector(FILE* ptrFile, E_Bool changeEndian, E_Bool formated,
                E_Int npts, char* varName, float*& fieldf, double*& fieldd)
{ 
  fieldf = NULL; fieldd = NULL;

  // recupere dataName, dataType
  char dataType[256];
  fscanf(ptrFile, " %s %s\n", varName, dataType);
  printf("varName=%s, dataType=%s\n", varName, dataType);

  if (strcmp(dataType, "float") == 0)
  {
    fieldf = new float [3*npts];
    if (formated == false)
    {
        fread(fieldf, sizeof(float), 3*npts, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%f %f %f", &fieldf[3*i], &fieldf[3*i+1], &fieldf[3*i+2]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < 3*npts; i++) fieldf[i] = FBE(fieldf[i]);
    }
  }
  else if (strcmp(dataType, "double") == 0)
  {
    fieldd = new double [3*npts];
    if (formated == false)
    {
        fread(fieldd, sizeof(double), 3*npts, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < npts; i++) fscanf(ptrFile, "%lf %lf %lf", &fieldd[3*i], &fieldd[3*i+1], &fieldd[3*i+2]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < 3*npts; i++) fieldd[i] = DBE(fieldd[i]);
    }
  }
}

//========================================================================
void readField(FILE* ptrFile, E_Bool changeEndian, E_Bool formated,
               E_Int& npts, E_Int& nfields, char* varName, float*& fieldf, double*& fieldd)
{ 
  fieldf = NULL; fieldd = NULL;

  // recupere dataName, numArrays
  char dataName[256]; E_Int numArrays;
  fscanf(ptrFile, " %s " SF_D_ "\n", dataName, &numArrays);
  printf("dataName=%s, numArrays=" SF_D_ "\n", dataName, numArrays);
  
  // Lit uniquement la premiere variable
  E_Int numComponents; E_Int numTuples; char dataType[256];
  fscanf(ptrFile, "%s " SF_D2_ " %s\n", varName, &numComponents, &numTuples, dataType);
  npts = numTuples;
  nfields = numComponents;
  E_Int size = npts*nfields;
  printf("varName=%s, dataType=%s, npts=" SF_D_ "\n", varName, dataType, npts);

  if (strcmp(dataType, "float") == 0)
  {
    fieldf = new float [size];
    if (formated == false)
    {
        fread(fieldf, sizeof(float), size, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < size; i++) fscanf(ptrFile, "%f", &fieldf[i]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < size; i++) fieldf[i] = FBE(fieldf[i]);
    }
  }
  else if (strcmp(dataType, "double") == 0)
  {
    fieldd = new double [size];
    if (formated == false)
    {
        fread(fieldd, sizeof(double), size, ptrFile);
        fgetc(ptrFile); // avoid \n
    }
    else
    {
        for (E_Int i = 0; i < size; i++) fscanf(ptrFile, "%lf", &fieldd[i]);
    }
    if (changeEndian)
    {
        for (E_Int i = 0; i < size; i++) fieldd[i] = DBE(fieldd[i]);
    }
  }
}

//=============================================================================
/* binvtkread
*/
//=============================================================================
E_Int K_IO::GenIO::binvtkread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: vtkread: cannot open file %s.\n", file);
    return 1;
  }

  E_Int i; char c;

  // Lecture de la version (jusqu'a newline)
  char buf[257];
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  printf("version: %s\n", buf);
  buf[22] = '\0';
  if (strcmp(buf, "# vtk DataFile Version") != 0) { fclose(ptrFile); return 1; }
  
  // Lecture du header (jusqu'a newline)
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i] = '\0';
  printf("header: %s\n", buf);
  
  // Type de format BINARY or ASCII
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i-1] = '\0';
  printf("FORMAT: %s\n", buf);
  E_Bool formated = true;
  E_Bool changeEndian = false;

  if (strcmp(buf, "BINARY") == 0) { formated = false; changeEndian = true; }
  printf("isformated=%d, changeEndian=%d\n", formated, changeEndian);
  
  // Type de Data (STRUCTURED_POINTS, STRUCTURED_GRID, RECTILINEAR_GRID, POLYDATA, UNSTRUCTURED_GRID)
  char dataSetType[256];
  i = 0; c = '1';
  while (c != '\n') { c = fgetc(ptrFile); buf[i] = c; i++; }
  buf[i-1] = '\0';
  E_Int l = strlen(buf);
  for (E_Int i = 0; i <= l-8; i++) dataSetType[i] = buf[i+8];
  dataSetType[i] = '\0';
  printf("DATA SET TYPE: %s\n", dataSetType);

  E_Int npts; float* ptsf; double* ptsd;
  int* cells; int* cellTypes; E_Int ncells; E_Int size;
    
  if (strcmp(dataSetType, "POLYDATA") == 0)
  {
    // Read POINTS
    readPoints(ptrFile, changeEndian, formated, npts, ptsf, ptsd);
    
    // READ ngon cells (only polygons for now)
    readPolygons(ptrFile, changeEndian, formated, ncells, size, cells);
  }
  else if (strcmp(dataSetType, "UNSTRUCTURED_GRID") == 0)
  {
    // Read POINTS
    readPoints(ptrFile, changeEndian, formated, npts, ptsf, ptsd);
    
    // READ CELLS
    readCells(ptrFile, changeEndian, formated, ncells, size, cells);
    readCellTypes(ptrFile, changeEndian, formated, ncells, cellTypes);
  }
  else 
  {
    printf("Warning: this kind of dataSet (%s) is not implemented.\n", dataSetType);
  }

  // READ POINT_DATA
  E_Int ret = readLine(ptrFile, buf);
  printf("keyword (POINT_DATA)=%s\n", buf);
  
  vector<float*> fieldsf;
  vector<double*> fieldsd;
  vector<E_Int> type; // 0: scalar, 1: vector
  char varName[256]; 
  varString = new char [1200];
  strcpy(varString, "x,y,z");
    
  // READ SCALARS, VECTORS, FIELD if ANY
  ret = fscanf(ptrFile, "%s", buf);  
  while (ret != EOF)
  {
    printf("keyword=%s\n", buf);
    if (strcmp(buf, "SCALARS") == 0)
    {
        float* ff; double* fd;
        readScalar(ptrFile, changeEndian, formated, npts, varName, ff, fd);
        strcat(varString, ",");
        strcat(varString, varName);
        if (ff != NULL) { fieldsf.push_back(ff); fieldsd.push_back(NULL); }
        else { fieldsd.push_back(fd); fieldsf.push_back(NULL); }
        type.push_back(1);
    }
    else if (strcmp(buf, "VECTORS") == 0)
    {
        float* ff; double *fd;
        readVector(ptrFile, changeEndian, formated, npts, varName, ff, fd);
        strcat(varString, ",");
        strcat(varString, varName);
        strcat(varString, "X,");
        strcat(varString, varName);
        strcat(varString, "Y,");
        strcat(varString, varName);
        strcat(varString, "Z");
        if (ff != NULL) { fieldsf.push_back(ff); fieldsd.push_back(NULL); }
        else { fieldsd.push_back(fd); fieldsf.push_back(NULL); }
        type.push_back(3);
    }
    else if (strcmp(buf, "FIELD") == 0)
    {
        float* ff; double* fd; E_Int nfields;
        readField(ptrFile, changeEndian, formated, npts, nfields, varName, ff, fd);
        if (nfields == 1)
        {
            strcat(varString, ",");
            strcat(varString, varName);
        }
        else if (nfields == 2)
        {
            for (E_Int k = 0; k < nfields; k++)
            strcat(varString, ",");
            strcat(varString, varName);
            strcat(varString, "X,");
            strcat(varString, varName);
            strcat(varString, "Y");
        }
        else if (nfields == 3)
        {
            strcat(varString, ",");
            strcat(varString, varName);
            strcat(varString, "X,");
            strcat(varString, varName);
            strcat(varString, "Y,");
            strcat(varString, varName);
            strcat(varString, "Z");
        }
        if (ff != NULL) { fieldsf.push_back(ff); fieldsd.push_back(NULL); }
        else { fieldsd.push_back(fd); fieldsf.push_back(NULL); }
        type.push_back(nfields);
    }
    ret = fscanf(ptrFile, "%s", buf);
  }

  // remise en array
  E_Int nvars = 0;
  for (size_t i = 0; i < type.size(); i++)
  {
    nvars += type[i];
  }
  
  // Concatenate arrays (fields supposed to be POINT_DATA)
  FldArrayF* f = new FldArrayF(npts, 3+nvars);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);

  if (ptsf != NULL)
  { 
    for (E_Int i = 0; i < npts; i++)
    {
        fx[i] = ptsf[3*i]; // cast
        fy[i] = ptsf[3*i+1];
        fz[i] = ptsf[3*i+2];
    }
    delete [] ptsf;
  }
  else
  {
    for (E_Int i = 0; i < npts; i++)
    {
        fx[i] = ptsd[3*i];
        fy[i] = ptsd[3*i+1];
        fz[i] = ptsd[3*i+2];
    }
    delete [] ptsd;
  }
    
  E_Int np = 4;
  for (size_t n = 0; n < type.size(); n++)
  {
    if (type[n] == 1)
    {
       E_Float* fp = f->begin(np);
       if (fieldsf[n] != NULL)
        for (E_Int i = 0; i < npts; i++) fp[i] = fieldsf[n][i];
       else for (E_Int i = 0; i < npts; i++) fp[i] = fieldsd[n][i];
       np += 1;
    }
    else if (type[n] == 2)
    {
       E_Float* fx = f->begin(np);
       E_Float* fy = f->begin(np+1);
    
       if (fieldsf[n] != NULL)
       {
        for (E_Int i = 0; i < npts; i++) fx[i] = fieldsf[n][2*i];
        for (E_Int i = 0; i < npts; i++) fy[i] = fieldsf[n][2*i+1];
       }
       else
       {
        for (E_Int i = 0; i < npts; i++) fx[i] = fieldsd[n][2*i];
        for (E_Int i = 0; i < npts; i++) fy[i] = fieldsd[n][2*i+1];
       }  
       np += 2; 
    }
    else if (type[n] == 3)
    {
       E_Float* fx = f->begin(np);
       E_Float* fy = f->begin(np+1);
       E_Float* fz = f->begin(np+2);
    
       if (fieldsf[n] != NULL)
       {
        for (E_Int i = 0; i < npts; i++) fx[i] = fieldsf[n][3*i];
        for (E_Int i = 0; i < npts; i++) fy[i] = fieldsf[n][3*i+1];
        for (E_Int i = 0; i < npts; i++) fz[i] = fieldsf[n][3*i+2];
       }
       else
       {
        for (E_Int i = 0; i < npts; i++) fx[i] = fieldsd[n][3*i];
        for (E_Int i = 0; i < npts; i++) fy[i] = fieldsd[n][3*i+1];
        for (E_Int i = 0; i < npts; i++) fz[i] = fieldsd[n][3*i+2];
       }  
       np += 3; 
    }
   }
    
   for (size_t n = 0; n < fieldsf.size(); n++) delete [] fieldsf[n];
   for (size_t n = 0; n < fieldsd.size(); n++) delete [] fieldsd[n];
    
   // Export de la connectivite
   if (strcmp(dataSetType, "UNSTRUCTURED_GRID") == 0)
   {
    // tri suivant les elements types
    E_Int nNODE = 0; // 1
    E_Int nBAR = 0; // 3
    E_Int nTRI = 0; // 5
    E_Int nQUAD = 0; // 9
    E_Int nHEXA = 0; // 12
    E_Int nTETRA = 0; // 10
    E_Int nPENTA = 0; // 13
    E_Int nPYRA = 0; // 14
    for (E_Int i = 0; i < ncells; i++)
    {
      // type de la cellule
      switch (cellTypes[i])
      {
        case 1:
          nNODE++; break;
        case 3:
          nBAR++; break;
        case 5: 
          nTRI++; break;
        case 9:
          nQUAD++; break;
        case 10:
          nTETRA++; break;
        case 12:
          nHEXA++; break;
        case 13:
          nPENTA++; break;
        case 14:
          nPYRA++; break;
        default: ;
      }
    }

    // Allocate
    printf("Elements: BAR=" SF_D_ " TRI=" SF_D_ " QUAD=" SF_D_ " HEXA=" SF_D_ " TETRA=" SF_D_ " PENTA=" SF_D_ " PYRA=" SF_D_ "\n", nBAR, nTRI, nQUAD, nHEXA, nTETRA, nPENTA, nPYRA);
    
    E_Int* cnBAR=NULL; E_Int* cnTRI=NULL; E_Int* cnQUAD=NULL;
    E_Int* cnTETRA=NULL; E_Int* cnPYRA=NULL; E_Int* cnPENTA=NULL;
    E_Int* cnHEXA=NULL;
    FldArrayI* cn;
    
    if (nBAR > 0) 
    { 
      cn = new FldArrayI(nBAR, 2); cnBAR = cn->begin(); 
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(1);
    }
    if (nTRI > 0)
    { 
      cn = new FldArrayI(nTRI, 3); cnTRI = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(2); 
    }
    if (nQUAD > 0) 
    { 
      cn = new FldArrayI(nQUAD, 4); cnQUAD = cn->begin(); 
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(3);
    }
    if (nHEXA > 0) 
    { 
      cn = new FldArrayI(nHEXA, 8); cnHEXA = cn->begin(); 
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(7);
    }
    if (nTETRA > 0) 
    { 
      cn = new FldArrayI(nTETRA, 4); cnTETRA = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(4); 
    }
    if (nPENTA > 0) 
    { 
      cn = new FldArrayI(nPENTA, 6); cnPENTA = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(6);
    }
    if (nPYRA > 0) 
    { 
      cn = new FldArrayI(nPYRA, 5); cnPYRA = cn->begin();
      unstructField.push_back(f); connect.push_back(cn); eltType.push_back(5);
    }

    int* ptr = cells;
    E_Int pBAR = 0; E_Int pPENTA = 0;
    E_Int pHEXA = 0; E_Int pTRI = 0; E_Int pQUAD = 0;
    E_Int pPYRA = 0; E_Int pTETRA = 0;
    for (E_Int i = 0; i < ncells; i++)
    {
      // type de la cellule
      switch (cellTypes[i])
      {
        case 3:
          cnBAR[pBAR] = ptr[1]+1;
          cnBAR[pBAR+nBAR] = ptr[2]+1;
          pBAR++;
          break;
        case 5:
          cnTRI[pTRI] = ptr[1]+1;
          cnTRI[pTRI+nTRI] = ptr[2]+1;
          cnTRI[pTRI+2*nTRI] = ptr[3]+1;
          pTRI++;
          break;
        case 9:
          cnQUAD[pQUAD] = ptr[1]+1;
          cnQUAD[pQUAD+nQUAD] = ptr[2]+1;
          cnQUAD[pQUAD+2*nQUAD] = ptr[3]+1;
          cnQUAD[pQUAD+3*nQUAD] = ptr[4]+1;
          pQUAD++;
          break;
        case 10:
          cnTETRA[pTETRA] = ptr[1]+1;
          cnTETRA[pTETRA+nTETRA] = ptr[2]+1;
          cnTETRA[pTETRA+2*nTETRA] = ptr[3]+1;
          cnTETRA[pTETRA+3*nTETRA] = ptr[4]+1;
          pTETRA++;
          break;
        case 12:
          cnHEXA[pHEXA] = ptr[1]+1;
          cnHEXA[pHEXA+nHEXA] = ptr[2]+1;
          cnHEXA[pHEXA+2*nHEXA] = ptr[3]+1;
          cnHEXA[pHEXA+3*nHEXA] = ptr[4]+1;
          cnHEXA[pHEXA+4*nHEXA] = ptr[5]+1;
          cnHEXA[pHEXA+5*nHEXA] = ptr[6]+1;
          cnHEXA[pHEXA+6*nHEXA] = ptr[7]+1;
          cnHEXA[pHEXA+7*nHEXA] = ptr[8]+1;
          pHEXA++;
          break;
        case 13:
          cnPENTA[pPENTA] = ptr[1]+1;
          cnPENTA[pPENTA+nPENTA] = ptr[2]+1;
          cnPENTA[pPENTA+2*nPENTA] = ptr[3]+1;
          cnPENTA[pPENTA+3*nPENTA] = ptr[4]+1;
          cnPENTA[pPENTA+4*nPENTA] = ptr[5]+1;
          cnPENTA[pPENTA+5*nPENTA] = ptr[6]+1;
          pPENTA++;
          break;
        case 14:
          cnPYRA[pPYRA] = ptr[1]+1;
          cnPYRA[pPYRA+nPYRA] = ptr[2]+1;
          cnPYRA[pPYRA+2*nPYRA] = ptr[3]+1;
          cnPYRA[pPYRA+3*nPYRA] = ptr[4]+1;
          cnPYRA[pPYRA+4*nPYRA] = ptr[5]+1;
          pPYRA++;
          break;
        default: ;
      }
      ptr = ptr+(ptr[0]+1);
    }
    
    delete [] cells;
    delete [] cellTypes;
   }
   else if (strcmp(dataSetType, "POLYDATA") == 0)
   {
       // export as a surface NGON
       E_Int sizeFN = 0; E_Int nfaces = 0;
       E_Int nelts = 0; E_Int sizeEF = 0;
       int* ptr = cells;
       E_Int n;
       for (E_Int i = 0; i < ncells; i++)
       {
           n = ptr[0];
           nfaces += n;
           sizeFN += 3*n;
           nelts += 1;
           sizeEF += n+1;
           ptr += n+1;
       }

       E_Int sizeTot = 4+sizeFN+sizeEF;
       FldArrayI* cn = new FldArrayI(sizeTot, 1); E_Int* cnp = cn->begin();
       ptr = cells;
       cnp[0] = nfaces;
       cnp[1] = sizeFN;
       cnp += 2;
       E_Int* cn2 = cnp+sizeFN;
       cn2[0] = nelts;
       cn2[1] = sizeEF;
       cn2 += 2;
       E_Int count = 1;
       for (E_Int i = 0; i < ncells; i++)
       {
           n = ptr[0];
           for (E_Int k = 0; k < n-1; k++) 
           {
            cnp[0] = 2;
            cnp[1] = ptr[k+1]+1;
            cnp[2] = ptr[k+2]+1;
            cnp += 3;
           }
           cnp[0] = 2;
           cnp[1] = ptr[n]+1; // loop
           cnp[2] = ptr[1]+1;
           cnp += 3;
           
           cn2[0] = n;
           for (E_Int k = 0; k < n; k++) 
           {
            cn2[k+1] = count+k;
           }
           count += n;
           cn2 += n+1;

           ptr += n+1;
       }
       unstructField.push_back(f); connect.push_back(cn); eltType.push_back(8);
       delete [] cells;
   }

  fclose(ptrFile);

  // Cree le nom de zone
  for (size_t i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }
  
  return 0;
}

//=============================================================================
// Only write NODE, BAR, TRI, QUAD, TETRA, HEXA, PYRA, PENTA.
//=============================================================================
E_Int K_IO::GenIO::binvtkwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // NODE, BAR, TRI, QUADS, TETRA, HEXA , PENTA, PYRA, supported
    if (eltTypes[zone][0] < 8) nvalidZones++;
    else
      printf("Warning: binvtkwrite: zone " SF_D_ " not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: binvtkwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  printf("Warning: binvtkwrite: not implemented.");
  return 1;
}

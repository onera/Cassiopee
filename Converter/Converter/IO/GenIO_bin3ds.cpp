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

// Binary 3DS (3D Studio) file support

# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include <stdio.h>
# include <string.h>

using namespace K_FLD; 
using namespace std;

//=============================================================================
/* 
   3dsread
*/
//=============================================================================
E_Int K_IO::GenIO::f3dsread( 
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
    printf("Warning: 3dsread: cannot open file %s.\n", file);
    return 1;
  }

  int i; // Index variable
  int obj = 0; // Nbre d'objets

  float tf;
  unsigned short ti;
  FldArrayF* f;
  E_Float *f1, *f2, *f3;
  FldArrayI* cp;
  E_Int *cp1, *cp2, *cp3;
  
  unsigned short l_chunk_id; // Chunk identifier
  unsigned int l_chunk_length; // Chunk length
  unsigned char l_char; // Char variable
  unsigned short l_qty; // Number of elements in each chunk

  unsigned short l_face_flags; // Flag that stores some face information

  // version
  int is3DS = 0;
  unsigned char version;
  KFSEEK(ptrFile, 28L, SEEK_SET);
  fread(&version, 1, 1, ptrFile);
  if (version < 3)
  { printf("Warning: 3dsread: this version is not really supported.\n"); }
  
  KFSEEK(ptrFile, 0, SEEK_END);
  E_LONG filelength = KFTELL(ptrFile);
  KFSEEK(ptrFile, 0, SEEK_SET);

  while (KFTELL(ptrFile) < filelength) // Loop to scan the whole file
  {
    fread (&l_chunk_id, 2, 1, ptrFile); // Read the chunk header 
    //printf("ChunkID: %x\n", l_chunk_id);
    fread (&l_chunk_length, 4, 1, ptrFile); // Read the length of the chunk
    //printf("ChunkLength: %d\n", l_chunk_length);
    if (l_chunk_length == 0) break;

    switch (l_chunk_id)
    {
      //----------------- MAIN3DS -----------------
      // Description: Main chunk, contains all the other chunks
      // Chunk ID: 4d4d
      // Chunk Length: 0 + sub chunks
      //-------------------------------------------
      case 0x4d4d:
        is3DS = 1;
        break;

        //----------------- EDIT3DS -----------------
        // Description: 3D Editor chunk, objects layout info
        // Chunk ID: 3d3d (hex)
        // Chunk Length: 0 + sub chunks
        //-------------------------------------------
      case 0x3d3d:
        break;

        //--------------- EDIT_OBJECT ---------------
        // Description: Object block, info for each object
        // Chunk ID: 4000 (hex)
        // Chunk Length: len(object name) + sub chunks
        //-------------------------------------------
      case 0x4000:
      {
        i = 0;
        obj++; 
        char* name = new char[128];
        //printf(" objet no: %i\n", obj);
        
        do
        {
          fread (&l_char, 1, 1, ptrFile);
          name[i] = l_char;
          i++;
        } while(l_char != '\0' && i < 127);
        name[i] = '\0';
        zoneNames.push_back(name);
        //printf(" objet: %s\n", name);
      }
      break;

        //--------------- OBJ_TRIMESH ---------------
        // Description: Triangular mesh, contains chunks for 3d mesh info
        // Chunk ID: 4100 (hex)
        // Chunk Length: 0 + sub chunks
        //-------------------------------------------
      case 0x4100:
        break;

        //--------------- TRI_VERTEX ---------------
        // Description: Vertices list
        // Chunk ID: 4110 (hex)
        // Chunk Length: 1 x unsigned short (number of vertices)
        // + 3 x float (vertex coordinates) x (number of vertices)
        // + sub chunks
        //-------------------------------------------
      case 0x4110:
        fread (&l_qty, sizeof(unsigned short), 1, ptrFile);
        f = new FldArrayF(l_qty, 3);
        unstructField.push_back(f);
        f1 = f->begin(1); f2 = f->begin(2); f3 = f->begin(3);
        //printf(" Number of vertices: %d\n", l_qty);
        for (i = 0; i < l_qty; i++)
        {
          fread (&tf, sizeof(float), 1, ptrFile); f1[i] = tf;
          fread (&tf, sizeof(float), 1, ptrFile); f2[i] = tf;
          fread (&tf, sizeof(float), 1, ptrFile); f3[i] = tf;
          //printf("Vertices %f %f %f\n", (*f)(i,1),(*f)(i,2), (*f)(i,3));
        }
        break;

        //--------------- TRI_FACEL1 ----------------
        // Description: Polygons (faces) list
        // Chunk ID: 4120 (hex)
        // Chunk Length: 1 x unsigned short (number of polygons)
        // + 3 x unsigned short (polygon points) x (number of polygons)
        // + sub chunks
        //-------------------------------------------
      case 0x4120:
        fread (&l_qty, sizeof(unsigned short), 1, ptrFile);
        //printf(" Number of polygons: %d\n", l_qty);
        eltType.push_back(2);
        cp = new FldArrayI(l_qty, 3);
        connect.push_back(cp);
        cp1 = cp->begin(1); cp2 = cp->begin(2); cp3 = cp->begin(3); 
        for (i = 0; i < l_qty; i++)
        {
          fread (&ti, sizeof(unsigned short), 1, ptrFile); cp1[i] = ti+1;
          fread (&ti, sizeof(unsigned short), 1, ptrFile); cp2[i] = ti+1;
          fread (&ti, sizeof(unsigned short), 1, ptrFile); cp3[i] = ti+1;
          fread (&l_face_flags, sizeof(unsigned short), 1, ptrFile);
          //printf("Face flags: %x\n", l_face_flags);
          //printf("%d %d %d\n", (*cp)(i,1), (*cp)(i,2), (*cp)(i,3));
        }
        break;

        //------------- TRI_MAPPINGCOORS ------------
        // Description: Vertices list
        // Chunk ID: 4140 (hex)
        // Chunk Length: 1 x unsigned short (number of mapping points)
        // + 2 x float (mapping coordinates) x (number of mapping points)
        // + sub chunks
        //-------------------------------------------
      case 0x4140:
        fread (&l_qty, sizeof (unsigned short), 1, ptrFile);
        for (i = 0; i < l_qty; i++)
        {
          fread (&tf, sizeof (float), 1, ptrFile);
          //printf("Mapping list u: %f\n",p_object->mapcoord[i].u);
          fread (&tf, sizeof (float), 1, ptrFile);
          //printf("Mapping list v: %f\n",p_object->mapcoord[i].v);
        }
        break;

        //----------- Skip unknown chunks ------------
        // We need to skip all the chunks that currently we don't use
        // We use the chunk length information to set the file pointer
        // to the same level next chunk
        //-------------------------------------------
      default:
        KFSEEK(ptrFile, l_chunk_length-6, SEEK_CUR);
    }
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);
  if (is3DS == 1) return 0;
  else return 1;
}

//=============================================================================
// Only write triangle arrays. Others are discarded.
//=============================================================================
E_Int K_IO::GenIO::f3dswrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  E_Int no = -1;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltType[zone] == 2) // triangles
    { nvalidZones++; if (no == -1) no = zone; }
    else
      printf("Warning: 3dswrite: zone %d not written (not a triangle zone).", zone);
  }

  if (nvalidZones == 0) return 1;
  
  // Zone must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: 3dswrite: zone do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  // Open file
  FILE* ptrFile = fopen(file, "wb");
  if (ptrFile == NULL) 
  {
    printf("3dswrite: I can't open file %s.\n", file);
    return 1;
  }
  
  // version + entete
  unsigned short chunk_id;
  unsigned int chunk_length;
  unsigned int til;
  unsigned short ti;

  // Size of main chunk (zones)
  E_Int size = 0;
  E_Int l4110, l4120, l4100, l4000;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltType[zone] == 2) // triangles
    {
      FldArrayF& coord = *unstructField[zone];
      E_Int nv = coord.getSize();
      FldArrayI& cn = *connect[zone];
      E_Int ne = cn.getSize();
      E_Int sizeVertex = 2 + 3*nv*4;
      E_Int sizeFace = 2 + 2*4*ne;
      l4110 = 6 + sizeVertex;
      l4120 = 6 + sizeFace;
      l4100 = 6 + l4110 + l4120;
      l4000 = 6 + strlen(zoneNames[zone]) + 1 + l4100;
      size = size + l4000;
    }
  }

  // MAIN3DS
  chunk_id = 0x4d4d;
  chunk_length = size + 10 + 6 + 6;
  fwrite(&chunk_id, 2, 1, ptrFile);
  fwrite(&chunk_length, 4, 1, ptrFile);

  // VERSION
  chunk_id = 0x0002;
  chunk_length = 10;
  fwrite(&chunk_id, 2, 1, ptrFile);
  fwrite(&chunk_length, 4, 1, ptrFile);
  til = 3; // version
  fwrite(&til, 4, 1, ptrFile);

  // EDIT3DS
  chunk_id = 0x3d3d;
  chunk_length = 6 + size;
  fwrite(&chunk_id, 2, 1, ptrFile);
  fwrite(&chunk_length, 4, 1, ptrFile);

  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltType[zone] == 2) // triangles
    {
      FldArrayF& coord = *unstructField[zone];
      E_Float* coordx = coord.begin(posx);
      E_Float* coordy = coord.begin(posy);
      E_Float* coordz = coord.begin(posz);
      E_Int nv = coord.getSize();
      FldArrayI& cn = *connect[zone];
      E_Int ne = cn.getSize();
      E_Int sizeVertex = 2 + 3*nv*4;
      E_Int sizeFace = 2 + 2*4*ne;
      l4110 = 6 + sizeVertex;
      l4120 = 6 + sizeFace;
      l4100 = 6 + l4110 + l4120;
      l4000 = 6 + strlen(zoneNames[zone]) + 1 + l4100;

      // EDIT object
      chunk_id = 0x4000;
      chunk_length = l4000;
      fwrite(&chunk_id, 2, 1, ptrFile);
      fwrite(&chunk_length, 4, 1, ptrFile);
      fwrite(zoneNames[zone], 1, strlen(zoneNames[zone])+1, ptrFile);

      // OBJECT TRIMESH
      chunk_id = 0x4100;
      chunk_length = l4100;
      fwrite(&chunk_id, 2, 1, ptrFile);
      fwrite(&chunk_length, 4, 1, ptrFile);

      // TRI_VERTEX      
      chunk_id = 0x4110;
      chunk_length = l4110;
      fwrite(&chunk_id, 2, 1, ptrFile);
      fwrite(&chunk_length, 4, 1, ptrFile);
      ti = nv;
      fwrite (&ti, sizeof(unsigned short), 1, ptrFile); // nv
      for (E_Int i = 0; i < nv; i++)
      {
        float x = coordx[i];
        float y = coordy[i];
        float z = coordz[i];
        
        fwrite(&x, 4, 1, ptrFile); 
        fwrite(&y, 4, 1, ptrFile);
        fwrite(&z, 4, 1, ptrFile);
      }

      // TRI_FACE
      chunk_id = 0x4120;
      chunk_length = l4120;
      fwrite(&chunk_id, 2, 1, ptrFile);
      fwrite(&chunk_length, 4, 1, ptrFile);
      ti = ne;
      fwrite (&ti, sizeof(unsigned short), 1, ptrFile); // ne
      for (E_Int i = 0; i < ne; i++)
      {
        ti = cn(i,1)-1;
        fwrite (&ti, sizeof(unsigned short), 1, ptrFile);
        ti = cn(i,2)-1;
        fwrite (&ti, sizeof(unsigned short), 1, ptrFile);
        ti = cn(i,3)-1;
        fwrite (&ti, sizeof(unsigned short), 1, ptrFile);
        ti = 0; // flags
        fwrite (&ti, sizeof(unsigned short), 1, ptrFile);
      }
    }
  }
  fclose(ptrFile);
  return 0;
}

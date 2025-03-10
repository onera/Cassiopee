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

// Binary gltf 2.0 (OpenGL) file support

# include "GenIO.h"
# include "Array/Array.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>
# define CGLTF_IMPLEMENTATION
# include "cgltf.h"

using namespace K_FLD;
using namespace std;

static void readAccessor(std::vector<float>& data, const cgltf_accessor* accessor)
{
  size_t components = cgltf_num_components(accessor->type);
  data.resize(accessor->count * components);
  cgltf_accessor_unpack_floats(accessor, &data[0], data.size());
}


// parse gltf data
void parseMeshesGltf(cgltf_data* data, std::vector<FldArrayF*>& unstructField, 
  std::vector<FldArrayI*>& connect, std::vector<E_Int>& eltType)
{
  cgltf_float lm[4*4];
  E_Int blockId = 0;
  // Pour tous les maillages
  for (size_t ni = 0; ni < data->nodes_count; ++ni)
  {
    cgltf_node& node = data->nodes[ni];

    if (!node.mesh) continue; // pas de maillage dans ce noeud

    const cgltf_mesh& mesh = *node.mesh;
    //int mesh_id = int(&mesh - data->meshes);
    int coordCreated = false;
    int texCoordFound = false;

    // transformation
    lm[ 0] = 1.; lm[ 1] = 0.; lm[ 2] = 0.; lm[ 3] = 0.;
    lm[ 4] = 0.; lm[ 5] = 1.; lm[ 6] = 0.; lm[ 7] = 0.;
    lm[ 8] = 0.; lm[ 9] = 0.; lm[10] = 1.; lm[11] = 0.;
    lm[12] = 0.; lm[13] = 0.; lm[14] = 0.; lm[15] = 1.;
    cgltf_node_transform_world(&node, lm);

    //printf("bloc: " SF_D_ ", mpc: %lld\n", blockId, mesh.primitives_count);
    
    // Pour tous ce qu'il y a dans mesh
    for (size_t pi = 0; pi < mesh.primitives_count; ++pi)
    {
      const cgltf_primitive& primitive = mesh.primitives[pi];
      if (primitive.type == cgltf_primitive_type_triangles)
      {
        //printf("bloc: " SF_D_ ", primitives are triangles\n", blockId);
        //size_t count = primitive.attributes ? primitive.attributes[0].data->count : 0;
        //for (size_t i = 0; i < count; ++i) printf("%zu\n", i);
      }
      if (primitive.type == cgltf_primitive_type_points)
      {
        printf("Warning: bingltf: bloc: " SF_D_ ", primitives are points\n", blockId);
      }
      
      if (primitive.type == cgltf_primitive_type_triangle_strip)
      {
        printf("Warning: bingltf: bloc: " SF_D_ ", primitives are triangle strips\n", blockId);
      }
      if (primitive.indices && primitive.type == cgltf_primitive_type_triangles)
      {
        //printf("bloc: " SF_D_ ", primitive are indices : %lld\n", blockId, primitive.indices->count);
        //printf("is this connect? %lld %lld\n", primitive.indices->count, primitive.indices->count/3);
        //for (size_t i = 0; i < primitive.indices->count; ++i)
        //  printf("%lld\n", cgltf_accessor_read_index(primitive.indices, i));
        
        E_Int nd = primitive.indices->count/3;
        eltType.push_back(2);
        FldArrayI* cp = new FldArrayI(nd, 3);
        connect.push_back(cp);
        E_Int* c1 = cp->begin(1);
        E_Int* c2 = cp->begin(2);
        E_Int* c3 = cp->begin(3);
        
        for (E_Int i = 0; i < nd; i++)
        {
          c1[i] = cgltf_accessor_read_index(primitive.indices, 3*i)+1;
          c2[i] = cgltf_accessor_read_index(primitive.indices, 3*i+1)+1;
          c3[i] = cgltf_accessor_read_index(primitive.indices, 3*i+2)+1;
        }
        //for (size_t i = 0; i < primitive.indices->count/3; ++i)
        //  printf("%lld\n", cgltf_accessor_read_index(primitive.indices, i));
      }
      
      // attributs
      for (size_t ai = 0; ai < primitive.attributes_count; ++ai)
      {
        const cgltf_attribute& attr = primitive.attributes[ai];

        if (attr.type == cgltf_attribute_type_invalid)
        {
          fprintf(stderr, "Warning: bingltf: ignoring unknown attribute %s in primitive %zu of mesh " SF_D_ "\n", attr.name, pi, blockId);
          continue;
        }

        // Lit cet attribut
        vector<float> data1;
        readAccessor(data1, attr.data);

        if (attr.type == cgltf_attribute_type_color && attr.data->type == cgltf_type_vec3)
        {
          //for (size_t i = 0; i < data.size(); ++i) data[i].f[3] = 1.0f;
        }

        if (attr.type == cgltf_attribute_type_position)
        {
          //printf("bloc " SF_D_ ": trouve coords\n", blockId);
          // numerotation i*3+j (i = pts, j = 0,1,2 composante)
          //for (E_Int i = 0; i < data.size(); i++) 
          //  printf("%f %f %f\n", data[i*3+0], data[i*3+1], data[i*3+2]);
          E_Int npts = data1.size()/3;
          FldArrayF* f;
          E_Float x,y,z;
          if (coordCreated) f = unstructField[blockId];
          else { f = new FldArrayF(npts, 8); unstructField.push_back(f); f->setAllValuesAtNull(); coordCreated = true; }
          E_Float* fx = f->begin(1);
          E_Float* fy = f->begin(2);
          E_Float* fz = f->begin(3);
          //for (E_Int i = 0; i < data.size(); i++) 
          //  printf("%f %f %f\n", data[i*3+0], data[i*3+1], data[i*3+2]);
          for (E_Int i = 0; i < npts; i++)
          {
            //printf("%f %f %f\n", data[i*3+0], data[i*3+1], data[i*3+2]);
            x = data1[i*3+0];
            y = data1[i*3+1];
            z = data1[i*3+2];

            //fx[i] = lm[0]*x + lm[1]*y + lm[2]*z + lm[3];
            //fy[i] = lm[5]*x + lm[5]*y + lm[6]*z + lm[7];
            //fz[i] = lm[8]*x + lm[9]*y + lm[10]*z + lm[11];
            fx[i] = lm[0]*x + lm[4]*y + lm[8]*z + lm[12];
            fy[i] = lm[1]*x + lm[5]*y + lm[9]*z + lm[13];
            fz[i] = lm[2]*x + lm[6]*y + lm[10]*z + lm[14];
          }
        }

        if (attr.type == cgltf_attribute_type_texcoord)
        {
          //printf("bloc " SF_D_ " : trouve texcoord\n", blockId);
          if (texCoordFound == true) continue; // to avoid multiple texcoord, the firts one is selected
          E_Int npts = data1.size()/2;
          FldArrayF* f;
          if (coordCreated) f = unstructField[blockId];
          else { f = new FldArrayF(npts, 8); unstructField.push_back(f); f->setAllValuesAtNull(); coordCreated = true; }
          E_Float* u = f->begin(4);
          E_Float* v = f->begin(5);
          for (E_Int i = 0; i < npts; i++)
          {
            //printf("%f %f\n", data[i*2+0], data[i*2+1]);
            u[i] = data1[i*2+0];
            v[i] = 1.-data1[i*2+1];
          }
          texCoordFound = true;
        }

        if (attr.type == cgltf_attribute_type_normal)
        {
          //printf("bloc " SF_D_ " : trouve normals\n", blockId);
          E_Int npts = data1.size()/3;
          FldArrayF* f;
          if (coordCreated) f = unstructField[blockId];
          else { f = new FldArrayF(npts, 8); unstructField.push_back(f); f->setAllValuesAtNull(); coordCreated = true; }
          E_Float* nx = f->begin(6);
          E_Float* ny = f->begin(7);
          E_Float* nz = f->begin(8);
        
          for (E_Int i = 0; i < npts; i++)
          {
            //printf("%f %f\n", data[i*3+0], data[i*3+1]);
            nx[i] = data1[i*3+0];
            ny[i] = data1[i*3+1];
            nz[i] = data1[i*3+2];
          }
        }

      }
   
      // morph target
      for (size_t ti = 0; ti < primitive.targets_count; ++ti)
      {
        const cgltf_morph_target& target = primitive.targets[ti];

        for (size_t ai = 0; ai < target.attributes_count; ++ai)
        {
          const cgltf_attribute& attr = target.attributes[ai];

          if (attr.type == cgltf_attribute_type_invalid)
          {
            fprintf(stderr, "Warning: bingltf: ignoring unknown attribute %s in morph target %zu of primitive %zu of mesh " SF_D_ "\n", attr.name, ti, pi, blockId);
            continue;
          }

          vector<float> data1;
          readAccessor(data1, attr.data);
        }
      }
      blockId++;
    }
  }
}

//=============================================================================
/* 
   gltfread 
*/
//=============================================================================
E_Int K_IO::GenIO::gltfread( 
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
    printf("Warning: gltfread: cannot open file %s.\n", file);
    return 1;
  }
  
  cgltf_data* data = 0;
  cgltf_options options = {};
  cgltf_result result = cgltf_parse_file(&options, file, &data);
  result = (result == cgltf_result_success) ? cgltf_load_buffers(&options, data, file) : result;
  result = (result == cgltf_result_success) ? cgltf_validate(data) : result;

  parseMeshesGltf(data, unstructField, connect, eltType);
  cgltf_free(data);
  
  // Cree les noms de zones
  for (E_Int i = 0; i < (E_Int)unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, i);
    zoneNames.push_back(zoneName);
  }
  varString = new char [30];
  strcpy(varString, "x,y,z,_u_,_v_,_nx_,_ny_,_nz_");
  
  return 0;
}

//=============================================================================
// Only write ONE triangle array. Others are discarded.
//=============================================================================
E_Int K_IO::GenIO::gltfwrite(
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
  E_Int no = -1;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    if (eltTypes[zone][0] == 2) // triangles
    { nvalidZones++; if (no == -1) no = zone; }
    else
      printf("Warning: gltfwrite: zone " SF_D_ " not written (not a triangle zone).", zone);
  } 

  if (nvalidZones == 0) return 1;
  if (nvalidZones > 1) printf("Warning: gltfwrite: only first zone will be written.");
  
  // Zone must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: gltfwrite: zone do not have coordinates. Not written.");
    return 1;
  }
  posx++; posy++; posz++;

  printf("Warning: gltfwrite: not implemented. Nothing written.\n");
  return 1;
}

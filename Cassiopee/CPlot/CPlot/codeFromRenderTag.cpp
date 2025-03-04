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

//=============================================================================
// transforme une color string en R,G,B
// color peut etre "White", ... ou en HEXA "#FFFFFF"
//=============================================================================
void Data::colorString2RGB(char* color, float& colorR, float& colorG, float& colorB)
{
  // Hard coded colors
  if (K_STRING::cmp(color, "White") == 0)
  { colorR = 1.; colorG = 1.; colorB = 1.; }
  else if (K_STRING::cmp(color, "Black") == 0)
  { colorR = 0.; colorG = 0.; colorB = 0.; }
  else if (K_STRING::cmp(color, "Grey") == 0)
  { colorR = 0.69; colorG = 0.69; colorB = 0.69; }
  else if (K_STRING::cmp(color, "Blue") == 0)
  { colorR = 73./255.; colorG = 86./255.; colorB = 243./255.; }
  else if (K_STRING::cmp(color, "Red") == 0)
  { colorR = 1.; colorG = 28./255.; colorB = 28./255.; }
  else if (K_STRING::cmp(color, "Green") == 0)
  { colorR = 7./255.; colorG = 180./255.; colorB = 3./255.; }
  else if (K_STRING::cmp(color, "Yellow") == 0)
  { colorR = 1.; colorG = 250./255.; colorB = 36./255.; }
  else if (K_STRING::cmp(color, "Orange") == 0)
  { colorR = 249./255.; colorG = 112./255.; colorB = 6./255.; }
  else if (K_STRING::cmp(color, "Magenta") == 0)
  { colorR = 230./255.; colorG = 0.0; colorB = 143./255.; }
  else if (K_STRING::cmp(color, "Brown") == 0)
  { colorR = 0.588; colorG = 0.294; colorB = 0.; }
  else if (color[0] == '#')  // par code hexa #aabbcc
  {
    char code[3]; unsigned int val;
    code[0] = color[1]; code[1] = color[2]; code[2] = '\0';
    sscanf(code, "%x", &val);
    colorR = val/255.;
    code[0] = color[3]; code[1] = color[4]; code[2] = '\0';
    sscanf(code, "%x", &val);
    colorG = val/255.;
    code[0] = color[5]; code[1] = color[6]; code[2] = '\0';
    sscanf(code, "%x", &val);
    colorB = val/255.;
  }
}

//=============================================================================
// IN: tag: chaine du tag
// OUT: parametres mis a jour a partir de la chaine tag
//=============================================================================
void Data::codeFromRenderTag(Zone& z, char* tag, 
                             float& colorR, float& colorG, float& colorB,
                             E_Int& material, double& blending, E_Int& meshOverlay,
                             float& shaderParam1, float& shaderParam2)
{
  colorR = -1.; colorG = -1.; colorB = -1.; material = -1; blending = -1.;
  meshOverlay = 0; shaderParam1 = 1.; shaderParam2 = 1.;

  if (tag == NULL) return;
  E_Int c = tag[0];
  char color[256];
  char mat[256];
  char temp[256];
  E_Int isoColor = 0;
  E_Int l;

  // Get color
  E_Int i = 0;
  while (c != '\0' && c != ':')
  {
    color[i] = c; i++; c = tag[i];
  }
  color[i] = '\0';
  if (K_STRING::cmp(color, "Iso") == 0)
  {
    isoColor = 1; l = 0;
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      color[l] = c; i++; l++; c = tag[i];
    }
    color[l] = '\0';
  }
  if (K_STRING::cmp(color, "centers") == 0 || K_STRING::cmp(color, "nodes") == 0)
  { 
    l = 0;
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      color[l] = c; i++; l++; c = tag[i];
    }
    color[l] = '\0';
  }

  // Get Material
  E_Int j = 0;
  if (c == '\0') // material missing
    strcpy(mat, "None");
  else
  {
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      mat[j] = c; j++; i++; c = tag[i];
    }
    mat[j] = '\0';
  }

  // Get blending
  j = 0;
  {
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      temp[j] = c; j++; i++; c = tag[i];
    }
    temp[j] = '\0';
  }
  if (K_STRING::cmp(temp, "None") != 0) blending = atof(temp);
  
  // Get meshOverlay
  j = 0;
  {
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      temp[j] = c; j++; i++; c = tag[i];
    }
    temp[j] = '\0';
  }
  if (K_STRING::cmp(temp, "None") != 0) meshOverlay = atoi(temp);

  // Get shader parameters
  j = 0;
  {
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      temp[j] = c; j++; i++; c = tag[i];
    }
    temp[j] = '\0';
  }
  if (K_STRING::cmp(temp, "None") != 0) shaderParam1 = atof(temp);

  j = 0;
  {
    i++; c = tag[i];
    while (c != '\0' && c != ':')
    {
      temp[j] = c; j++; i++; c = tag[i];
    }
    temp[j] = '\0';
  }
  if (K_STRING::cmp(temp, "None") != 0) shaderParam2 = atof(temp);
  //printf("shader parameters: %f %f\n", shaderParam1, shaderParam2);

  // Analyse couleur
  if (isoColor == 0) 
  {
    colorString2RGB(color, colorR, colorG, colorB);
  }
  else
  { // iso color
    colorR = 1.; colorG = 0.; colorB = 0.;
    for (E_Int n = 0; n < z.nfield; n++)
    {
      if (K_STRING::cmp(z.varnames[n], color) == 0) { colorR = -n-2; break;}
    }
  }

  // Analyse material
  material = 0; // None or Solid
  if (K_STRING::cmp(mat, "Glass") == 0) material = 1;
  else if (K_STRING::cmp(mat, "Chrome") == 0) material = 2;
  else if (K_STRING::cmp(mat, "Metal") == 0) material = 3;
  else if (K_STRING::cmp(mat, "Wood") == 0) material = 4;
  else if (K_STRING::cmp(mat, "Marble") == 0) material = 5;
  else if (K_STRING::cmp(mat, "Smoke") == 0) material = 6;
  else if (K_STRING::cmp(mat, "XRay") == 0) material = 7;
  else if (K_STRING::cmp(mat, "Granite") == 0) material = 8;
  else if (K_STRING::cmp(mat, "Sphere") == 0) material = 9;
  else if (K_STRING::cmp(mat, "Brick") == 0) material = 10;
  else if (K_STRING::cmp(mat, "Cloud") == 0) material = 11;
  else if (K_STRING::cmp(mat, "Gooch") == 0) material = 12;
  else if (K_STRING::cmp(mat, "Flat") == 0) material = 13;
  else if (K_STRING::cmp(mat, "Texmat") == 0) material = 14;  
}

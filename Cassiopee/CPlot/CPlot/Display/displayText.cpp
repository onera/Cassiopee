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

// Routines pour afficher des textes sur l'ecran

#include "Data.h"

// Choix de la font bitmap.
// 0: Fonts GLUT
// 1: Fonts OpenGLText (non fonctionnel, a finaliser)
#define FONTMETHOD 0

//=============================================================================
// Choix pour la fonction de render de texte :
// renderBitmapString1 : avec les polices de glut
// renderBitmapString2 : avec le texte en texture
//=============================================================================
void Data::renderBitmapString(float x, float y, float z,
                              E_Int fontSize, char *string,
                              float colorR, float colorG, float colorB, float colorA,
                              float nx, float ny, float nz,
                              float r)
{
#if FONTMETHOD == 0
  renderBitmapString1(x, y, z, fontSize, string, colorR, colorG, colorB, colorA, nx, ny, nz);
#else
  renderBitmapString2(x, y, z, fontSize, string, colorR, colorG, colorB, colorA, nx, ny, nz);
#endif
}
//=============================================================================
// Choix pour la fonction mesure de la largeur du texte :
// textWidth1 : avec les polices de glut
// textWidth2 : avec le texte en texture
//=============================================================================
E_Int Data::textWidth(E_Int fontSize, char* string)
{
#if FONTMETHOD == 0
  return textWidth1(fontSize, string);
#else
  return textWidth2(fontSize, string);
#endif
}

//=============================================================================
// Retourne la font GLUT
//=============================================================================
void* Data::getGlutFont(E_Int fontSize)
{

  if (fontSize == 10) return GLUT_BITMAP_HELVETICA_10;
  else if (fontSize == 12) return GLUT_BITMAP_HELVETICA_12;
  else if (fontSize == 18) return GLUT_BITMAP_HELVETICA_18;
  else return GLUT_BITMAP_HELVETICA_12;
}
//=============================================================================
// Retourne la largeur en pixels de la chaine
//=============================================================================
E_Int Data::textWidth1(E_Int fontSize, char* string)
{
  void* font = getGlutFont(fontSize);
  E_Int l = 0;
  char* c;
  for (c = string; *c != '\0'; c++) 
  {
    l += glutBitmapWidth(font, *c);
  }
  return l;
}
//=============================================================================
// Retourne la hauteur en pixels de la police
//=============================================================================
E_Int Data::textHeight(E_Int fontSize)
{
  return fontSize;
}

//=============================================================================
// Render a bitmap string at a given position
// IN: xyz: position of string
// IN: font: type de font
// IN: myString: chaine a afficher
// IN: nx,ny,nz: vecteur vertical du plan d'ecriture (norme)
// IN: r: ratio de distance
// Don't deal with \n
//=============================================================================
void Data::renderBitmapString1(float x, float y, float z,
                               E_Int fontSize, char *myString,
                               float colorR, float colorG, float colorB, float colorA,
                               float nx, float ny, float nz,
                               float r) 
{
  void* font = getGlutFont(fontSize);
  glColor4f(colorR, colorG, colorB, colorA);
  glRasterPos3f(x, y, z);
  E_Int i = 0;
  while (myString[i] != '\0')
  {
    glutBitmapCharacter(font, myString[i]); i++;
  }
}

//==============================================================================
// Render a bitmap string at a given position with Shadow
// Don't deal with \n, deal with color specified with @1, @2, ...
// IN: xyz: position of string
// IN: fontSize: taille de la font
// IN: myString: chaine
// IN: fgColor: couleur ecriture
// IN: shColor: couleur shadow
// IN: offt: vecteur horizontal du plan d'ecriture
// IN: offn: vecteur vertical du plan d'ecriture
// IN: r: ratio de distance (avec shadow)
//==============================================================================
void Data::renderStringWithShadow(
  float x, float y, float z,
  E_Int fontSize, char *myString,
  float fgColorR, float fgColorG, float fgColorB, float fgColorA,
  float shColorR, float shColorG, float shColorB, float shColorA,
  double offtx, double offty, double offtz,
  double offnx, double offny, double offnz,
  double r) 
{
  // local colors
  float lColorR; float lColorG; float lColorB; 
  lColorR = fgColorR; lColorG = fgColorG; lColorB = fgColorB;  

  double shnx, shny, shnz, shtx, shty, shtz;
  
  // decalage pour le cas avec depth test (bandeaux)
  double zoffset = 0.;
  if (offtx == 1 && offty == 0 && offtz == 0) zoffset=0.01;
  
  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  shtx = offtx*r; shty = offty*r; shtz = offtz*r;
  shnx = offnx*r; shny = offny*r; shnz = offnz*r;
  
  // segmente la string en mots
  char msg[1024];
  E_Int i = 0; E_Int j = 0; E_Int w = 0;
  while (myString[j] != '\0')
  {
    if (myString[j] == '@')
    {
      // Ecrit
      msg[i] = '\0';
      renderBitmapString(x+w+shnx+shtx, y+shny+shty, z+shnz+shtz-zoffset, 
                         fontSize, msg, shColorR, shColorG, shColorB, shColorA);
      renderBitmapString(x+w, y, z, fontSize, msg, lColorR, lColorG, lColorB, fgColorA);
      w += textWidth(fontSize, msg); i = 0;

      // Change la couleur
      j++;
      if (myString[j] == '0') // default
      { lColorR = fgColorR; lColorG = fgColorG; lColorB = fgColorB; }
      if (myString[j] == '2') // light green
      { lColorR = 252./255.; lColorG = 113./255.; lColorB = 1./255.; }
      else if (myString[j] == '1') // light red 
      { lColorR = 155./255.; lColorG = 251./255.; lColorB = 60./255.; }
    
      j++;
    }

    /*
    if (myString[j] == ' ')
    {
      msg[i] = '\0';
      glColor4f(lColorR, lColorG, lColorB, fgColorA);
      renderBitmapString(x+w, y, z, fontSize, msg);
      glColor4f(shColorR, shColorG, shColorB, shColorA);
      renderBitmapString(x+w+shnx+shtx, y+shny+shty, z+shnz+shtz, fontSize, msg);
      renderBitmapString(x+w+shtx, y+shty, z+shtz, fontSize, msg);
      renderBitmapString(x+w+shnx, y+shny, z+shnz, fontSize, msg);
      w += textWidth(fontSize, msg); i = 0;
      lColorR = fgColorR; lColorG = fgColorG; lColorB = fgColorB; 
    }
    */
    msg[i] = myString[j]; i++; j++;
  }
  msg[i] = '\0';
  renderBitmapString(x+w+shnx+shtx, y+shny+shty, z+shnz+shtz-zoffset, 
                    fontSize, msg, shColorR, shColorG, shColorB, shColorA);
  renderBitmapString(x+w, y, z, fontSize, msg, lColorR, lColorG, lColorB, fgColorA);
  //glDisable(GL_BLEND);
}

//=============================================================================
// Put back the previous perspective
//=============================================================================
void Data::resetPerspectiveProjection() 
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}

//=============================================================================
// Set orthographic projection
//=============================================================================
void Data::setOrthographicProjection() 
{
  // Switch to projection mode
  glMatrixMode(GL_PROJECTION);
  // Save previous matrix which contains the 
  // settings for the perspective projection
  glPushMatrix();
  // Reset matrix
  glLoadIdentity();
  // Set a 2D orthographic projection
  gluOrtho2D(0, _view.w, 0, _view.h);
  // Invert the y axis, down is positive
  glScalef(1, -1, 1);
  // Move the origin from the bottom left corner
  // to the upper left corner
  glTranslatef(0, -_view.h, 0);
  glMatrixMode(GL_MODELVIEW);
}

//=============================================================================
// Display the text at the upper corner (font1)
//=============================================================================
void Data::displayText(char* text)
{
  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix(); glLoadIdentity();

  // Render a string
  renderStringWithShadow(5, _font1Size, 0, _font1Size, text,
                         1., 1., 1., 1.,
                         0.1, 0.1, 0.1, 0.5);

  // Render a rectangle
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.0, 0.0, 1., 0.5);
  glBegin(GL_QUADS);
  glVertex3d(0, 0, 0);
  glVertex3d(_view.w, 0, 0);
  glVertex3d(_view.w, _font1Size+5, 0);
  glVertex3d(0, _font1Size+5, 0);
  glEnd();
  glDisable(GL_BLEND);

  // Put back the previous projection
  glPopMatrix();
  resetPerspectiveProjection(); 
}

//=============================================================================
// Display big text (font2)
//=============================================================================
void Data::displayBigText(E_Int posx, E_Int posy, char* text)
{
  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix(); glLoadIdentity();

  // Render a string
  renderBitmapString(posx, posy, 0, _font2Size, text, 1., 1., 1., 1.);

  // Put back the previous projection
  glPopMatrix();
  resetPerspectiveProjection(); 
}

//=============================================================================
// Display small text (font1)
//=============================================================================
void Data::displaySmallText(E_Int posx, E_Int posy, char* text)
{
  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();

  // Render a string
  renderBitmapString(posx, posy, 0, _font1Size, text, 1., 1., 1., 1.);
  glPopMatrix();

  // Put back the previous projection
  resetPerspectiveProjection(); 
}

//=============================================================================
// Display CPlot message header
//=============================================================================
void Data::printHeader()
{
  char msg[256]; char tmp[256];

  // CPlot
  strcpy(msg, "CPlot");
  strcat(msg, " | ");
  // Dim Mode
  if (ptrState->dim != 1)
  {
    // 2D or 3D mode
    if (ptrState->dim == 2)
    {
      strcat(msg, "2D Mode (");
      switch (ptrState->var2D)
      {
        case 0:
          strcat(msg, "x,y) | ");
          break;
        case 1:
          strcat(msg, "x,z) | ");
          break;
        case 2:
          strcat(msg, "y,z) | ");
          break;
      }
    }
    else
      strcat(msg, "3D Mode | ");
  }
  else
  {
    // 1D mode
    if (ptrState->ijk1D == 0) strcat(msg, "1Di Mode (");
    else if (ptrState->ijk1D == 1) strcat(msg, "1Dj Mode (");
    else strcat(msg, "1Dk Mode (");

    switch (ptrState->var1D)
    {
      case 0:
        strcat(msg, "x,f) | ");
        break;
      case 1:
        strcat(msg, "y,f) | ");
        break;
      case 2:
        strcat(msg, "z,f) | ");
        break;
      case 3:
        strcat(msg, "s,f) | ");
        break;
    }
  }
  // Display mode
  strcat(msg, "Display ");
  switch (ptrState->mode)
  {
    case MESH:
      strcat(msg, "mesh");
      break;

    case SOLID:
      strcat(msg, "solid");
      break;
      
    case RENDER:
      strcat(msg, "render");
      break;

    case SCALARFIELD:
      if (_numberOfZones > 0)
        strcat(msg, _zones[0]->varnames[ptrState->scalarField]);
      break;

    case VECTORFIELD:
      if (_numberOfZones > 0)
      {
        strcat(msg, _zones[0]->varnames[ptrState->vectorField1]);
        strcat(msg, ", ");
        strcat(msg, _zones[0]->varnames[ptrState->vectorField2]);
        strcat(msg, ", ");
        strcat(msg, _zones[0]->varnames[ptrState->vectorField3]);
      }
      break;
  }

  if (_niter > 0)
  {
    strcat(msg, " | It. "); sprintf(tmp, SF_D_, _niter);
    strcat(msg, tmp);
  }

  if (ptrState->selectedZone > 0 && ptrState->selectedZone <= _numberOfZones)
  {
    int nz = ptrState->selectedZone-1;

    strcat(msg, " | ");
    strcat(msg, _zones[nz]->zoneName);
    strcat(msg, " : ");
    if (nz < _numberOfStructZones)
    {
      sprintf(tmp, SF_D_, _szones[nz]->ni);
      strcat(msg, tmp);
      strcat(msg, "x");
      sprintf(tmp, SF_D_, _szones[nz]->nj);
      strcat(msg, tmp);
      strcat(msg, "x");
      sprintf(tmp, SF_D_, _szones[nz]->nk);
      strcat(msg, tmp);

      // Displayed planes
      strcat(msg, " |");
      if (ptrState->dim != 1)
      {
        if (_szones[nz]->iPlane >= 0)
        {
          strcat(msg, " i=");
          sprintf(tmp, SF_D_, _szones[nz]->iPlane+1);
          strcat(msg, tmp);
        }
        if (_szones[nz]->jPlane >= 0)
        {
          strcat(msg, "+j=");
          sprintf(tmp, SF_D_, _szones[nz]->jPlane+1);
          strcat(msg, tmp);
        }
        if (_szones[nz]->kPlane >= 0)
        {
          strcat(msg, "+k=");
          sprintf(tmp, SF_D_, _szones[nz]->kPlane+1);
          strcat(msg, tmp);
        }
      }
      else
      {
        switch (ptrState->ijk1D)
        {
          case 0: // i variant
            if (_szones[nz]->jPlane >= 0)
            {
              strcat(msg, " j=");
              sprintf(tmp, SF_D_, _szones[nz]->jPlane+1);
              strcat(msg, tmp);
            }
            strcat(msg, " k=");
            sprintf(tmp, SF_D_, _szones[nz]->kLine+1);
            strcat(msg, tmp);
            break;
          
          case 1: // j variant
            strcat(msg, " i=");
            sprintf(tmp, SF_D_, _szones[nz]->iLine+1);
            strcat(msg, tmp);
            if (_szones[nz]->kPlane >= 0)
            {
              strcat(msg, " k=");
              sprintf(tmp, SF_D_, _szones[nz]->kPlane+1);
              strcat(msg, tmp);
            }
            break;
            
          case 2: // k variant
            if (_szones[nz]->iPlane >= 0)
            {
              strcat(msg, " i=");
              sprintf(tmp, SF_D_, _szones[nz]->iPlane+1);
              strcat(msg, tmp);
            }
            strcat(msg, " j=");
            sprintf(tmp, SF_D_, _szones[nz]->jLine+1);
            strcat(msg, tmp);
            break;
        }
      }
    }
    else // non structure
    {
      sprintf(tmp, SF_D_ " pts, " SF_D_ " elts", _uzones[nz-_numberOfStructZones]->np, 
              _uzones[nz-_numberOfStructZones]->ne);
      strcat(msg, tmp);
    }
  }

  displayText(msg);
}

//=============================================================================
// Display temporary message on line two (font2)
//=============================================================================
void Data::printTmpMessage(const char* text)
{
  char txt[512];
  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();
  E_Int length = strlen(text);
  E_Int c = 0; E_Int d = 0;
  E_Int l = 0;
  while (c < length)
  {
    // Build string
    d = 0;
    while (text[c] != '\0' && text[c] != '\n')
    {
      txt[d] = text[c]; c++; d++;
    }
    txt[d] = '\0'; c++;

    // Render string
    renderStringWithShadow(5, _font2Size+5+(l+1)*_font2Size, 0, _font2Size, txt,
                           1.0, 1.0, 1.0, 1.0,
                           0.1, 0.1, 0.1, 0.5);
    l++;
  }
  glPopMatrix();
  resetPerspectiveProjection();
}

//=============================================================================
// Display info window
//=============================================================================
void Data::displayInfoWindow(char* text, E_Int l)
{
  char msg[1024]; char msg2[1024]; E_Int j; E_Int jl;
  // Swap to orthographic 2D projection view
  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();

  // Render a string with newline
  E_Int ls = l-20;
  E_Int i1 = 0; // on text
  E_Int i2 = 0; // on msg
  char c = text[i1];
  E_Int w = 0;
  while (c != '\0')
  {
    if (c == '\n')
    {
      msg[i2] = '\0';
      // suppress @
      jl = 0;
      for (j = 0; j <= i2; j++)
      {
        if (msg[j] == '@') j++;
        else { msg2[jl] = msg[j]; jl++; }
      }
      
      w = MAX(w, textWidth(_font1Size, msg2));
      renderStringWithShadow(5, l, 0, _font1Size, msg,
                             1.0, 1.0, 1.0, 1.0,
                             0.1, 0.1, 0.1, 0.5);
      i2 = 0; l += 17;
    }
    else
    {
      msg[i2] = text[i1];
      i2++;
    }
    i1++;
    c = text[i1];
  }
  msg[i2] = '\0';
  // suppress @
  jl = 0;
  for (j = 0; j <= i2; j++)
  {
    if (msg[j] == '@') j++;
    else { msg2[jl] = msg[j]; jl++; }
  }
  w = MAX(w, textWidth(_font1Size, msg2));
  renderStringWithShadow(5, l, 0, _font1Size, msg,
                         1.0, 1.0, 1.0, 1.0,
                         0.1, 0.1, 0.1, 0.5);
  l += 16;

  // Render rectangle
  w += 10;
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.0, 0.0, 1., 0.5);
  glBegin(GL_QUADS);
  glVertex3d(2, l-5, 0);
  glVertex3d(w, l-5, 0);
  glVertex3d(w, ls, 0);
  glVertex3d(2, ls, 0);
  glEnd();
  glDisable(GL_BLEND);
  glPopMatrix();

  // Put back the previous projection
  resetPerspectiveProjection(); 
}

//=============================================================================
// Display informations over the display
//=============================================================================
void Data::displayInfo()
{
  char msg[1024]; char temp[512]; char local[512];

  // Nombre de champs
  E_Int nv = 0;
  if (_numberOfZones > 0) nv = _zones[0]->nfield;

  // Nombre de lignes a afficher
  E_Int nlignes = 0;

  // Display info on position (cam) (if no selected point)
  sprintf(msg, " ");
  if (ptrState->selectedZone <= 0)
  {
    sprintf(msg, "Cam: %g, %g, %g\nEye: %g, %g, %g\nDir: %g, %g, %g",  
            _view.xcam, _view.ycam, _view.zcam,
            _view.xeye, _view.yeye, _view.zeye,
            _view.dirx, _view.diry, _view.dirz);
    nlignes += 3;
  }

  // Display info on active zone
  if (ptrState->selectedZone > 0)
  {
    Zone* z = _zones[ptrState->selectedZone-1];
    
    if (ptrState->selectedZone-1 < _numberOfStructZones)
    {
      StructZone* zs = (StructZone*)z;
      sprintf(temp,"@1%s@0 (STRUCT): " SF_D_ "x" SF_D_ "x" SF_D_ "",
              z->zoneName, zs->ni, zs->nj, zs->nk);
      strcat(msg, temp);
    }
    else
    {
      UnstructZone* zu = (UnstructZone*)z;
      E_Int ncon = ptrState->activePointL;
      E_Int eltType = zu->eltType[ncon];
      E_Int eltSize = zu->eltSize[ncon];
      E_Int* connect = zu->connect[ncon];
      switch (eltType)
      {
        case 0:
          sprintf(temp,"@1%s@0 (NODE): " SF_D_ " pts",
                  z->zoneName, zu->npts);
          break;
        case 1:
          if (zu->_is_high_order && eltSize == 3)
            sprintf(temp,"@1%s@0 (BAR_3): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (BAR): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 2:
          if (zu->_is_high_order && eltSize == 6)
            sprintf(temp,"@1%s@0 (TRI_6): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (TRI): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 3:
          if (zu->_is_high_order && eltSize == 8)
            sprintf(temp,"@1%s@0 (QUAD_8): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else if (zu->_is_high_order && eltSize == 9)
            sprintf(temp,"@1%s@0 (QUAD_9): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (QUAD): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 4:
          if (zu->_is_high_order && eltSize == 10)
            sprintf(temp,"@1%s@0 (TETRA_10): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (TETRA): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 5:
          if (zu->_is_high_order && eltSize == 15)
            sprintf(temp,"@1%s@0 (PENTA_15): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else if (zu->_is_high_order && eltSize == 18)
            sprintf(temp,"@1%s@0 (PENTA_18): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (PENTA): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 6:
          if (zu->_is_high_order && eltSize == 14)
            sprintf(temp,"@1%s@0 (PYRA_14): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (PYRA): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 7:
          if (zu->_is_high_order && eltSize == 20)
            sprintf(temp,"@1%s@0 (HEXA_20): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else if (zu->_is_high_order && eltSize == 27)
            sprintf(temp,"@1%s@0 (HEXA_27): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          else
            sprintf(temp,"@1%s@0 (HEXA): " SF_D_ " pts, " SF_D_ " elts",
                    z->zoneName, zu->npts, zu->ne);
          break;
        case 10:
          sprintf(temp,"@1%s@0 (NGON): " SF_D_ " pts, " SF_D_ " elts, " SF_D_ " faces",
                  z->zoneName, zu->npts, zu->ne, connect[0]);
          break;
        default:
          sprintf(temp,"@1%s@0 (UNKNOWN): " SF_D_ " pts, " SF_D_ " elts",
                  z->zoneName, zu->npts, zu->ne);
          break;
      }
      strcat(msg, temp);
    }
    nlignes++;
   
    if (ptrState->selectedZone-1 < _numberOfStructZones)
    {
      // structure
      if (ptrState->mode == MESH || ptrState->mode == SOLID)
        strcpy(local, "\nx=@1%g@0, y=@1%g@0, z=@1%g@0 (@1i=" SF_D_ "@0, @1j=" SF_D_ "@0, @1k=" SF_D_ "@0)");
      else strcpy(local, "\nx=%g, y=%g, z=%g (i=" SF_D_ ", j=" SF_D_ ",  k=" SF_D_ ")");
      sprintf(temp, local,
              ptrState->activePointX, ptrState->activePointY, ptrState->activePointZ,
              ptrState->activePointI, ptrState->activePointJ, abs(ptrState->activePointK));
    }
    else
    {
      // non structure
      // np: no du vertex
      // ne: no de l'element
      // ns: no du sommet correpondant a np dans ne
      // ou nf: no de la face (NGON)
      Zone* z = _zones[ptrState->selectedZone-1];
      UnstructZone* zu = (UnstructZone*)z;
      if (zu->eltType[0] != 10)
      {
        if (ptrState->mode == MESH || ptrState->mode == SOLID)
          strcpy(local, "\nx=@1%g@0, y=@1%g@0, z=@1%g@0 (@1np=" SF_D_ "@0, @1ne=" SF_D_ "@0, @1ns=" SF_D_ "@0)");
        else strcpy(local, "\nx=%g, y=%g, z=%g (np=" SF_D_ ", ne=" SF_D_ ", ns=" SF_D_ ")");
      }
      else
      {
        if (ptrState->mode == MESH || ptrState->mode == SOLID)
          strcpy(local, "\nx=%g, y=%g, z=%g (@1np=" SF_D_ "@0, @1ne=" SF_D_ "@0, @1nf=" SF_D_ "@0)");
        else strcpy(local, "\nx=%g, y=%g, z=%g (np=" SF_D_ ", ne=" SF_D_ ", nf=" SF_D_ ")");
      }
      
      sprintf(temp, local,  
              ptrState->activePointX, ptrState->activePointY, ptrState->activePointZ,
              ptrState->activePointI, ptrState->activePointJ, abs(ptrState->activePointK));
    }
    strcat(msg, temp);
    nlignes++;

    if (nv > 0)
    {
      strcat(msg, "\n");
      
      E_Int count = 0;
      
      if (nv < 10) // moins de 10 variables, on affiche tout
      {
        for (E_Int i = 0; i < nv; i++)
        {
          if (i == ptrState->scalarField && ptrState->mode != MESH && ptrState->mode != SOLID)
            sprintf(temp, "@1%s=%g@0",  
                    z->varnames[i], ptrState->activePointF[i]);
          else
            sprintf(temp, "%s=%g",  
                    z->varnames[i], ptrState->activePointF[i]);
          strcat(msg, temp);
          count += strlen(temp);
          if (i < nv-1) strcat(msg, ", ");
          if (count > 50 && i < nv-1)
          {nlignes++; strcat(msg, "\n"); count = 0;}
        }
        nlignes++;
      }
      else // on affiche le debut et autour de la variable courante
      {
        E_Int tp[10];
        if (ptrState->scalarField < 2)
        {
          tp[0] = 0; tp[1] = 1; tp[2] = 2; tp[3] = 3; tp[4] = 4; tp[5] = 5;
          tp[6] = -1; 
          tp[7] = nv-3; tp[8] = nv-2; tp[9] = nv-1;
        }
        else if (ptrState->scalarField < nv-2)
        {
          tp[0] = 0; tp[1] = 1; tp[2] = 2; 
          tp[3] = -1; 
          tp[4] = ptrState->scalarField-1; tp[5] = ptrState->scalarField; tp[6] = ptrState->scalarField+1;
          tp[7] = -1; 
          tp[8] = nv-2; tp[9] = nv-1;
        }
        else
        {
          tp[0] = 0; tp[1] = 1; tp[2] = 2; tp[3] = 3; tp[4] = 4; tp[5] = -1;
          tp[6] = nv-4; 
          tp[7] = nv-3; tp[8] = nv-2; tp[9] = nv-1;
        }
        for (E_Int i = 0; i < 10; i++)
        {
          E_Int tt = tp[i];
          if (tt != -1)
          {
            if (tp[i] == ptrState->scalarField && ptrState->mode != MESH && ptrState->mode != SOLID)
              sprintf(temp, "@1%s=%g@0",  
                      z->varnames[tt], ptrState->activePointF[tt]);
            else
              sprintf(temp, "%s=%g",  
                      z->varnames[tt], ptrState->activePointF[tt]);
            strcat(msg, temp);
          }
          else 
          {
            strcpy(temp, "..."); strcat(msg, temp);
          }
          
          count += strlen(temp);
          if (i < 9) strcat(msg, ", ");
          if (count > 50 && i < nv-1)
          {nlignes++; strcat(msg, "\n"); count = 0;}
        }
        nlignes++;
      }
    }
  }
  displayInfoWindow(msg, _view.h-15*nlignes-4);
}

//============================================================================
// Create openGLText si pas encore cree
// Il y a deux tailles de fontes _font1Size et _font2Size
// en 1280x720 : 12, 18
// en 1920x1080: 18, 27
// en 3840x2160: 36, 54
// a l'ecran, je crois que windows upscale. Il faut donc utiliser les 
// polices 12, 18
// Par contre, pour les sorties offscreen, il faut utiliser les polices
// correspondants a la resolution de l'export.
//============================================================================
OpenGLText* Data::getOpenGLText(E_Int fontSize)
{
  char fontName[256];
  OpenGLText* pt = NULL;
  if (fontSize == 12)
  { 
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/BRITANIC_12");
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/lucon_12");
    strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/consola_12");
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/shelley_12");
  }
  else if (fontSize == 18) 
  {
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/BRITANIC_18");
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/lucon_18");
    strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/consola_18");
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/shelley_18");
  }
  else 
  {
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/BRITANIC_12");
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/lucon_12");
    strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/consola_12");
    //strcpy(fontName,  "../../Apps/Modules/CPlot/CPlot/Fonts/shelley_12");
  }
  if (fontSize == _font1Size) pt = _oglText1;
  else if (fontSize == _font2Size) pt = _oglText2;
  else pt = _oglText3;
  
  if (pt != NULL) return pt;
  E_Int canvasWidth = _view.w; E_Int canvasHeight = _view.h;
  pt = new OpenGLText();
  
  if (fontSize == _font1Size) _oglText1 = pt;
  else if (fontSize == _font2Size) _oglText2 = pt;
  else _oglText3 = pt;

  if (!pt->init(fontName, canvasWidth, canvasHeight)) return NULL;
  return pt;
}

//=============================================================================
// Render a bitmap string at a given position
// IN: xyz: position of string
// IN: font: type de font
// IN: myString: chaine a afficher
// IN: nx,ny,nz: vecteur vertical du plan d'ecriture (norme)
// IN: r: ratio de distance
// Manque : gestion des couleurs, gestion de l'affichage dans l'espace physique
//=============================================================================
void Data::renderBitmapString2(float x, float y, float z,
                               E_Int fontSize, char *myString,
                               float colorR, float colorG, float colorB, float colorA,
                               float nx, float ny, float nz,
                               float r) 
{
  // Load font
  OpenGLText* pt = getOpenGLText(fontSize);
  OpenGLText& oglText = *pt;
  
  E_Int canvasWidth = _view.w; E_Int canvasHeight = _view.h;
  oglText.changeSize(canvasWidth, canvasHeight);
  oglText.changeCanvas(canvasWidth, canvasHeight);
  
  oglText.backupStates();

  // render string
  E_Int posX = x; E_Int posY = _view.h-1-y;
  oglText.beginString();
  //float bbStr[2];
  //oglText.stringSize(myString, bbStr);
  float color[4];
  color[0] = colorR; color[1] = colorG; color[2] = colorB; color[3] = colorA;  
  oglText.drawString(posX, posY, myString, 0, color);
  oglText.endString(); // will render the whole at once
  oglText.restoreStates();
}

//=============================================================================
// Retourne la largeur en pixels de la chaine
//=============================================================================
E_Int Data::textWidth2(E_Int fontSize, char* string)
{
  OpenGLText* pt = getOpenGLText(fontSize);
  float bbStr[2];
  pt->stringSize(string, bbStr);
  return bbStr[0];
}

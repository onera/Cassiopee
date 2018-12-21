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
#include "../Data.h"
#include <stdlib.h>
#include "gl2ps.h"

//=============================================================================
// Screen dump plugins
//=============================================================================

//=============================================================================
// Fait le rendu, dump et ecrit le fichier
//=============================================================================
void Data::exportFile()
{
  if ((ptrState==NULL) || (ptrState->_isExporting==1)) {
    pthread_mutex_lock(&ptrState->export_mutex);
    //if (d->ptrState->shootScreen == 1)
    pthread_cond_wait (&ptrState->unlocked_export, &ptrState->export_mutex); 
    pthread_mutex_unlock(&ptrState->export_mutex);    
  }
  int stateHeader, stateInfo, stateMenu, stateBB;
  stateHeader = ptrState->header;
  stateInfo = ptrState->info;
  stateMenu = ptrState->menu;
  stateBB = ptrState->bb;
  ptrState->header = 0;
  ptrState->info = 0;
  ptrState->menu = 0;
  ptrState->bb = 0;
  dumpWindow();
  ptrState->header = stateHeader;
  ptrState->info = stateInfo;
  ptrState->menu = stateMenu;
  ptrState->bb = stateBB;
  if (ptrState->continuousExport == 0) { ptrState->shootScreen = 0; }  
  ptrState->_mustExport = 1;
  pthread_cond_signal(&ptrState->unlocked_export); // signal end of export
  pthread_mutex_unlock(&ptrState->export_mutex);
}

//=============================================================================
// Display to an image using FBO
// Retourne un buffer contenant l'image RGB
// exportWidth et exportHeight doivent etre Pairs
// Si mode=1, ajoute depth au buffer
//=============================================================================
char* Data::export2Image(int exportWidth, int exportHeight) 
{
  //printf("mode = %d\n", ptrState->offscreen);

  // resolution
  GLuint fb, rb, db;
#ifdef __SHADERS__
  glGenFramebuffersEXT(1, &fb);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fb);
  // Create and attach a color buffer
  glGenRenderbuffersEXT(1, &rb);
  // We must bind color_rb before we call glRenderbufferStorageEXT
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rb);
  // The storage format is RGBA8
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA8, 
                           exportWidth, exportHeight);
  // Attach color buffer to FBO
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
                               GL_RENDERBUFFER_EXT, rb);
  
  glGenRenderbuffersEXT(1, &db);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, db);
  if (ptrState->offscreen > 2) 
  {
    //glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH24_STENCIL8, exportWidth, exportHeight);
    //glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT32, exportWidth, exportHeight);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT32F, exportWidth, exportHeight);
    // Attach depth buffer to FBO
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER_EXT, db);
  }
  else 
  {
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, exportWidth, exportHeight);
    // Attach depth buffer to FBO
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                                 GL_RENDERBUFFER_EXT, db);
  }

  int status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (status == GL_FRAMEBUFFER_COMPLETE_EXT)
  {
    // SUCCESS
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fb);
  }
  else 
  { // FAIL
    printf("Warning: CPlot: unable to get the correct frame buffer.\n");
    exportWidth = _view.w; exportHeight = _view.h;
    exportWidth = (exportWidth/2)*2; exportHeight = (exportHeight/2)*2;
  }
#endif

  int s = exportWidth*exportHeight;
  unsigned szDepth = 4*s; // Nombre d'octets par pixel si DEPTH32F et STENCIL 8 : 4 bytes pour detpth + 1 byte  spare + 1 byte stencil
  // _accumBuffer = 4*s pixels à afficher + 4*s donnée de profondeur :-)
  //if (_accumBuffer.size() == 0) std::vector<unsigned char>(8*s, unsigned char(0)).swap(_accumBuffer);
  char* buffer;
  buffer = (char*)malloc(s*3*sizeof(char));
  float* depth = (float*)malloc(szDepth);
  for (E_Int i = 0; i < s; i++) depth[i] = 0;

  // Switch viewport
  int viewWSav = _view.w; int viewHSav = _view.h;
  _view.w = exportWidth; _view.h = exportHeight;
  glViewport(0, 0, (GLsizei) exportWidth, (GLsizei) exportHeight);
  _view.ratio = (double)_view.w/(double)_view.h;

  if (ptrState->stereo == 0) display();
  else displayAnaglyph();

  // Back viewport
  _view.w = viewWSav; _view.h = viewHSav;
  _view.ratio = (double)_view.w/(double)_view.h;
  glViewport(0, 0, (GLsizei) _view.w, (GLsizei) _view.h);

  // Traduction en buffer
  if (ptrState->offscreen != 1)
  {
    if (ptrState->offscreen > 2)
    {
      if (ptrState->offscreenBuffer == NULL)
      {
        //float scale, bias;
        //glGetFloatv(GL_DEPTH_SCALE, &scale);
        //glGetFloatv(GL_DEPTH_BIAS, &bias);
        //printf("scale=%f bias=%f\n", scale, bias);
        
        glReadPixels(0, 0, exportWidth, exportHeight, GL_RGB, GL_UNSIGNED_BYTE, buffer);

        glReadPixels(0, 0, exportWidth, exportHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depth);
        //glReadPixels(0, 0, exportWidth, exportHeight, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT_24_8, depth);
        //glReadPixels(0, 0, exportWidth, exportHeight, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, depth);
        //glReadPixels(0, 0, exportWidth, exportHeight, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, depth);
        //glReadPixels(0, 0, exportWidth, exportHeight, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, depth);

        //glGetFloatv(GL_DEPTH_SCALE, &scale);
        //glGetFloatv(GL_DEPTH_BIAS, &bias);GL_UNSIGNED_INT_24_8
        //printf("scale=%f bias=%f\n", scale, bias);
        
        ptrState->offscreenBuffer = (char*)malloc(s*3*sizeof(char));
        memcpy(ptrState->offscreenBuffer, buffer, s*3*sizeof(char));
        ptrState->offscreenDepthBuffer = (float*)malloc(szDepth);
        /*GLfloat projMat[16];
        glGetFloatv(GL_PROJECTION_MATRIX, projMat);
        for ( int i = 0; i < 4; ++i ) {
          for ( int j = 0; j < 4; ++j )
            printf("%11.8f\t", projMat[4*j+i]);
          printf("\n");
        }
        */
        double zNear = _view.nearD; 
        double zFar  = _view.farD;
        //printf("Znear : %e, zFar : %e\n",zNear, zFar);
        for (int i = 0; i < exportHeight*exportWidth; i++)
        {
          double z_n = 2.*double(depth[i])-1.0;
          double z_e = 2.0*zNear*zFar/(zFar+zNear-z_n*(zFar-zNear));
          ptrState->offscreenDepthBuffer[i] = float(z_e);//0.5*(-A*depth[i]+B)/depth[i] + 0.5;
        }
        
        //memcpy(ptrState->offscreenDepthBuffer, depth, szDepth);
        /*printf("Profondeur obtenue pour le rendu 1:\n");
        for ( size_t i = 0; i < exportHeight; ++i ) {
          for ( size_t j = 0; j < exportWidth; ++j )
            printf("%9.7e",ptrState->offscreenDepthBuffer[i*exportWidth+j]);
          printf("\n");
        }*/
        
      }
      else
      {
                // Relit et stocke RGB + depth
        glReadPixels(0, 0, exportWidth, exportHeight, GL_RGB, GL_UNSIGNED_BYTE, buffer);
        glReadPixels(0, 0, exportWidth, exportHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depth);
        /*GLfloat projMat[16];
        glGetFloatv(GL_PROJECTION_MATRIX, projMat);
        for ( int i = 0; i < 4; ++i ) {
          for ( int j = 0; j < 4; ++j )
            printf("%11.8f", projMat[4*j+i]);
          printf("\n");
        }*/
        double zNear = _view.nearD; 
        double zFar  = _view.farD;
        //printf("Znear : %e, ZFar : %e\n",zNear, zFar);
        for (int i = 0; i < exportHeight*exportWidth; i++)
        {
          double z_n = 2.*double(depth[i])-1.0;
          double z_e = 2.0*zNear*zFar/(zFar+zNear-z_n*(zFar-zNear));
          depth[i] = float(z_e);//0.5*(-A*depth[i]+B)/depth[i] + 0.5;
        }
        /*printf("Profondeur obtenue pour le rendu 2:\n");
        for ( size_t i = 0; i < exportHeight; ++i ) {
          for ( size_t j = 0; j < exportWidth; ++j )
            printf("%9.7e",depth[i*exportWidth+j]);
          printf("\n");
        }*/
        //exit(0);

        char* offscreen = (char*)ptrState->offscreenBuffer;
        float* offscreenD = (float*)ptrState->offscreenDepthBuffer;
        for ( unsigned i = 0; i < exportHeight; ++i ) {
          for ( unsigned j = 0; j < exportWidth; ++j ) {
            unsigned ind = i*exportWidth+j;
            assert(ind<s);
            if (depth[ind] <offscreenD[ind]) {
              //printf("replace pixel %d,%d => %e < %e \n", i, j, depth[ind], offscreenD[ind]);
              offscreen[3*ind  ] = buffer[3*ind  ];
              offscreen[3*ind+1] = buffer[3*ind+1];
              offscreen[3*ind+2] = buffer[3*ind+2];
              offscreenD[ind]    = depth[ind];
            }
          }
        }
        memcpy(buffer, ptrState->offscreenBuffer, 3*s*sizeof(char));
      }
    }
    else
    {
      if (ptrState->offscreenBuffer != NULL)
          memcpy(buffer, ptrState->offscreenBuffer, 3*s*sizeof(char));
      else
          glReadPixels(0, 0, exportWidth, exportHeight, GL_RGB, GL_UNSIGNED_BYTE, buffer);
    }
  }
  else 
  { // mesa offscreen rendering
    char* buffRGBA = (char*)ptrState->offscreenBuffer;
    for (int i = 0; i < s; i++)
    { 
      buffer[3*i] = buffRGBA[4*i];
      buffer[3*i+1] = buffRGBA[4*i+1];
      buffer[3*i+2] = buffRGBA[4*i+2];
    }
  }

#ifdef __SHADERS__
  // Delete FBO
  glDeleteRenderbuffersEXT(1, &rb);
  glDeleteRenderbuffersEXT(1, &db);
  //Bind 0, which means render to back buffer, as a result, fb is unbound
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  glDeleteFramebuffersEXT(1, &fb);
#endif

  free(depth);
  return buffer;
}

//=============================================================================
/*
  Dump a window to a buffer suivant l'extension du plugins.
*/
//=============================================================================
void Data::dumpWindow() 
{
  char fileName[120];
  if (_pref.screenDump == NULL) return;

  // File name
  if (strcmp(ptrState->exportFile, "CPlot") == 0)
    sprintf(fileName, "%s/%s.%d.%s", ptrState->localPathName, 
            ptrState->exportFile, ptrState->exportNumber, 
            _pref.screenDump->extension);
  else strcpy(fileName, ptrState->exportFile);
  ptrState->exportNumber++;

  if (strcmp(_pref.screenDump->extension, "ps") == 0)
  {
    // Postscript
    int state = GL2PS_OVERFLOW, buffsize = 0;

    FILE* fp = fopen(fileName, "wb");
    while (state == GL2PS_OVERFLOW)
    {
      buffsize += 1024*1024;
      gl2psBeginPage("test", "CPlot", NULL, 
                     GL2PS_EPS, GL2PS_BSP_SORT, 
		     GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT, 
		     GL_RGBA, 0, NULL, 0, 0, 0,  buffsize, fp, fileName);
      display();
      state = gl2psEndPage();
    }
    printf("Wrote file %s.\n", fileName);
    fclose(fp);
  }
  else // Other formats
  {

    int exportWidth = _view.w;
    int exportHeight = _view.h;
    double r = _view.h * 1. / _view.w;
    if (ptrState->exportWidth != -1 && ptrState->exportHeight != -1)
    {
      exportWidth = ptrState->exportWidth; 
      exportHeight = int(exportWidth * r);
    }
    else if (ptrState->exportWidth != -1 && ptrState->exportHeight == -1)
    {
      exportWidth = ptrState->exportWidth;
      exportHeight = int(exportWidth * r);
    }
    else if (ptrState->exportWidth == -1 && ptrState->exportHeight != -1)
    {
      exportHeight = ptrState->exportHeight;
      exportWidth = int(exportHeight * (1./r));
    }
    exportWidth = (exportWidth/2)*2;
    exportHeight = (exportHeight/2)*2; // doit etre pair
    
    char* buffer; char* buffer2;
    int antialiasing = 1;
    if (ptrState->offscreen >= 3) antialiasing=0;

    if (antialiasing == 1)
    {
      // get image X2
      buffer = export2Image(2*exportWidth, 2*exportHeight);
    
      // blur image X2
      buffer2 = (char*)malloc(3*2*exportWidth*2*exportHeight*sizeof(char));
      int nitBlur;
      if (ptrState->mode == MESH) nitBlur = int(exportWidth/500.)+1;
      else nitBlur = 1;
      gaussianBlur(2*exportWidth, 2*exportHeight, buffer, buffer2, nitBlur, 0.1);

      // supersample X2
      free(buffer);
      buffer = (char*)malloc(3*exportWidth*exportHeight*sizeof(char));
      superSample(exportWidth, exportHeight, buffer2, buffer, 2);
      free(buffer2);
    }
    else
    {
      buffer = export2Image(exportWidth, exportHeight);
    }
    
    if (strcmp(_pref.screenDump->extension, "mpeg") == 0)
    {
      buffer2 = (char*)malloc(3*exportWidth*exportHeight*sizeof(char));
      gaussianBlur(exportWidth, exportHeight, buffer, buffer2, 1, 0.1);
      //sharpenImage(exportWidth, exportHeight, buffer, buffer2, 5.,
      //             2, 2);
      for (int i = 0; i < 3*exportWidth*exportHeight; i++) buffer[i] = buffer2[i];
      free(buffer2);
    }

    // Dump the buffer to a file
    if (ptrState->offscreen != 3)
      _pref.screenDump->f(this, fileName, buffer, exportWidth, exportHeight, 0);    
    free(buffer);
  }
}

//=============================================================================
void Data::finalizeExport()
{
  int exportWidth = _view.w;
  int exportHeight = _view.h;
  if (ptrState->exportWidth != -1) exportWidth = ptrState->exportWidth;
  if (ptrState->exportHeight != -1) exportHeight = ptrState->exportHeight;
  exportWidth = (exportWidth/2)*2;
  exportHeight = (exportHeight/2)*2; // doit etre pair
  char* buffer = NULL;
  // Dump the buffer to a file
  _pref.screenDump->f(this, ptrState->exportFile, buffer, 
                      exportWidth, exportHeight, 1);
}

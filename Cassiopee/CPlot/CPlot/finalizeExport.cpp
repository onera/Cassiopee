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
# include "cplot.h"
# include "Data.h"

PyObject* K_CPLOT::isDisplayRunning(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  long isRunning = long(d->_CDisplayIsLaunched);
  return Py_BuildValue("l", isRunning);
}
//=============================================================================
/* Arrete l'export continu, finalize les MPEG, fait une barriere pour
   attendre la fin de l'ecriture */
//=============================================================================
PyObject* K_CPLOT::finalizeExport(PyObject* self, PyObject* args)
{
  int finalizeType = 0;
  if (!PyArg_ParseTuple(args, "i", &finalizeType)) return NULL;

  Data* d = Data::getInstance();
  
  // Finalize pour osmesa (delete le context)
  if (finalizeType == 1 || finalizeType == 5 || finalizeType == 6 || finalizeType == 7)
  {
#ifdef __MESA__
    // We may sometimes need to delete the context (if new image if a new size)
    // But generally, this is not the case
    //free(d->ptrState->offscreenBuffer[d->ptrState->frameBuffer]);
    //d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = NULL;
    //OSMesaDestroyContext(*(OSMesaContext*)(d->ptrState->ctx));
    //d->ptrState->ctx = NULL;
    if (finalizeType == 6)
    {
      // in composite mode, this buffer must be deleted
      //free(d->ptrState->offscreenBuffer[d->ptrState->frameBuffer]);
      //d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = NULL;
      free(d->ptrState->offscreenBuffer[d->ptrState->frameBuffer+1]);
      d->ptrState->offscreenBuffer[d->ptrState->frameBuffer+1] = NULL;
      free(d->ptrState->offscreenDepthBuffer[d->ptrState->frameBuffer]);
    }
#endif
    return Py_BuildValue("l", KSUCCESS);
  }

  // Autre finalize
  if ((d->ptrState == NULL) || (d->ptrState->_mustExport == 0)) 
  {
    pthread_mutex_lock(&d->ptrState->export_mutex);
    //if (d->ptrState->shootScreen == 1)
    pthread_cond_wait (&d->ptrState->unlocked_export, &d->ptrState->export_mutex); 
    pthread_mutex_unlock(&d->ptrState->export_mutex);    
  }
  d->ptrState->_isExporting = 1;
  // Bloc en attendant la fin de l'ecriture
  //if (d->ptrState->continuousExport == 0)
  //{}
  
  // Finalize mpeg
  if (finalizeType == -1 && strcmp(d->_pref.screenDump->extension, "mpeg") == 0)
    d->finalizeExport(); // force l'ecriture finale du fichier
  
  d->ptrState->continuousExport = 0;
  d->ptrState->shootScreen = 0;
  d->ptrState->_mustExport = 0;
  d->ptrState->_isExporting = 0;

  if (finalizeType == 4)
  {
    // clear compositing buffers
    free(d->ptrState->offscreenBuffer[d->ptrState->frameBuffer]);
    d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = NULL;
    free(d->ptrState->offscreenDepthBuffer[d->ptrState->frameBuffer]);
  }
  pthread_cond_signal(&d->ptrState->unlocked_export); // signal end of export
  pthread_mutex_unlock(&d->ptrState->export_mutex);

  return Py_BuildValue("l", KSUCCESS);
}

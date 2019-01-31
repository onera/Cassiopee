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
clearDisplay();
glLoadIdentity();
E_Float dx = _view.xeye-_view.xcam;
E_Float dy = _view.yeye-_view.ycam;
E_Float dz = _view.zeye-_view.zcam;
E_Float dirx = _view.dirx;
E_Float diry = _view.diry;
E_Float dirz = _view.dirz;
E_Float nx = dy*dirz-dz*diry;
E_Float ny = dz*dirx-dx*dirz;
E_Float nz = dx*diry-dy*dirx;
E_Float norm = MAX(sqrt(dirx*dirx + diry*diry + dirz*dirz), 1.e-10);
dirx = dirx / norm; diry = diry / norm; dirz = dirz / norm;
norm = MAX(sqrt(nx*nx + ny*ny + nz*nz), 1.e-10);
nx = nx / norm; ny = ny / norm; nz = nz / norm;
norm = MAX(sqrt(dx*dx + dy*dy + dz*dz), 1.e-10);
dx = dx / norm; dy = dy / norm; dz = dz / norm;
nx = ptrState->lightOffsetX*nx-ptrState->lightOffsetY*dirx+dx; 
ny = ptrState->lightOffsetX*ny-ptrState->lightOffsetY*diry+dy; 
nz = ptrState->lightOffsetX*nz-ptrState->lightOffsetY*dirz+dz;
E_Float norm2 = MAX(sqrt(nx*nx+ny*ny+nz*nz), 1.e-10);
nx = nx / norm2; ny = ny / norm2 ; nz = nz / norm2;
gluLookAt(
  _view.xeye-1.3*norm*nx, 
  _view.yeye-1.3*norm*ny, 
  _view.zeye-1.3*norm*nz, 
  _view.xeye+norm*nx, 
  _view.yeye+norm*ny, 
  _view.zeye+norm*nz, 
  _view.dirx, _view.diry, _view.dirz);

// Save light matrices in matrix texture
glGetDoublev(GL_MODELVIEW_MATRIX, _lightModelView);
glGetDoublev(GL_PROJECTION_MATRIX, _lightProjection);
glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // render only depth
glEnable(GL_CULL_FACE);
glCullFace(GL_FRONT);
//glEnable(GL_POLYGON_OFFSET_FILL);
//glPolygonOffset(1.1, 4);
if (ptrState->dim != 1)
{
  displaySSolid(); displayUSolid();
}
if (_shadowMap != 0) glDeleteTextures(1, &_shadowMap);
glActiveTexture(GL_TEXTURE0);
glGenTextures(1, &_shadowMap);
glBindTexture(GL_TEXTURE_2D, _shadowMap);
glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 0, 0, 
                 _view.w, _view.h, 0);
glDisable(GL_CULL_FACE);
//glPolygonOffset(0., 0.);
//glDisable(GL_POLYGON_OFFSET_FILL);
glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE); // render all

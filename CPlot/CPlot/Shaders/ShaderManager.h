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
#ifndef _CPLOT_SHADERMANAGER_HPP_
#define _CPLOT_SHADERMANAGER_HPP_
#include <cassert>
#include <map>
#include <vector>

#include <iostream>
#include "Shader.h"
#include "TesselationShaderManager.hpp"

// prepare the shadow texture
#define SHADOWTEXTURE if (ptrState->shadow > 0) {          \
  glActiveTexture(GL_TEXTURE0);                         \
  glBindTexture(GL_TEXTURE_2D, _shadowMap);                             \
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);    \
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);    \
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);          \
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);          \
  glMatrixMode(GL_TEXTURE); glLoadIdentity(); glLoadMatrixd(_bias);     \
  glMultMatrixd(_lightProjection); glMultMatrixd(_lightModelView);      \
  glMatrixMode(GL_MODELVIEW); }

namespace CPlot
{
  class ShaderManager
  {
  public:
    ShaderManager();
    ~ShaderManager();
    
    /* init GLEW/ARB */
    int init();
    
    /* load les shaders */
    int load();

    /* Add a new vertex/fragment shader from a vertex and 
       fragment files.
       This shader will be compiled and linked
       This subroutine returns an indentificator for the new shader object 
    */
    Shader* addFromFile(const char* vertexFile, const char* fragmentFile);
    Shader* addFromFile(const char* geomFile, const char* vertexFile, const char* fragmentFile);

    bool eraseShader(Shader* obj);
    unsigned short numberOfShaders() const
    { return (unsigned short)_shaderList.size(); }
    unsigned short currentShader() const
    { return _currentActiveShader; }
    unsigned short shader_id( int idShadMat ) const
    {
      unsigned short id = idShadMat*(tesselationManager.numberOfShaders()+1) + tesselationManager.currentShader();
      if ( id >= (unsigned short)_shaderList.size() ) std::cerr << "Bogue, depassement de tableau de shader !!!!" << std::flush << std::endl;
      return id;
    }

    void set_tesselation( unsigned short idTes );
    void unset_tesselation();
    bool has_tesselation()const { return tesselationManager.is_activate(); }
    
    Shader* operator [] (unsigned short id)
    {
        assert(id > 0);
        assert(id < _shaderList.size());
        return _shaderList[id]; 
    }
    
    const Shader* operator [] (unsigned short id) const
    { 
        assert(id > 0);
        assert(id < _shaderList.size());
        return _shaderList[id]; 
    }
    
    /* Return the id of a shader in shader list */
    unsigned short getId(Shader* shad) const;
    
    /* Activate a shader from his identificator 
       Desactivate the previous shader if different.
    */
    void activate(unsigned short id);
    /* Activate a shader (if contained in Shadermanager)
       Desactivate the previous shader if different.
    */
    void activate(Shader* shad);
    /* Deactivate shader pipeline and returns to a
       default opengl pipeline ( sauf si tesselation actif... )
    */
    void deactivate();
    
  private:
    ShaderManager(const ShaderManager& shadMan);
    ShaderManager& operator = (const ShaderManager& shadMan);
    std::vector<Shader*> _shaderList;
    Shader* m_previous_shader;
    unsigned short _currentActiveShader;
    TesselationShaderManager tesselationManager;
  };
}
#endif

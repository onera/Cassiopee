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
#ifndef _CPLOT_SHADER_HPP_
#define _CPLOT_SHADER_HPP_
#include "GL/glew.h"
#include <list>
#include <string>
#include "ShaderObject.h"
#include "FragmentShader.h"

namespace CPlot
{
  class Shader
  {
  public:
    Shader();
    virtual ~Shader();

    void add(const ShaderObject& shad);
    
    GLhandleARB getProgramId() const { return _programId; }

    bool start(); // Activate the shader. If not linked, link this shader before
    std::string getLinkerLog() const; // Return error log of the linker.
    void end(); // Desactivate the shader. OpenGL calls with the standard pipeline.

    unsigned int getAttributeLocation(const char* attributeName) const;
    unsigned int getUniformLocation  (const char* uniformName  ) const;

    // Submitting Uniform Variables. 
    // To apply change in shader after changing uniform variables,
    // you must stop and restart your shader (Open GL requirement)
    void setUniform(const char* varname, float v0);
    void setUniform(unsigned int idVar, float v0);
 
    void setUniform(const char* varname, float v0, float v1);
    void setUniform(unsigned int idVar, float  v0, float v1);

    void setUniform(const char* varname, float v0, float v1, float v2);
    void setUniform(unsigned int idVar , float v0, float v1, float v2);

    void setUniform(const char* varname, float v0, float v1, float v2, float v3);
    void setUniform(unsigned int idVar , float v0, float v1, float v2, float v3);
    
    void setUniform(const char* varname, int v0);
    void setUniform(unsigned int idVar , int v0);

    void setUniform(const char* varname, int v0, int v1);
    void setUniform(unsigned int idVar , int v0, int v1);

    void setUniform(const char* varname, int v0, int v1, int v2);
    void setUniform(unsigned int idVar , int v0, int v1, int v2);

    void setUniform(const char* varname, int v0, int v1, int v2, int v3);
    void setUniform(unsigned int idVar , int v0, int v1, int v2, int v3);
    
    // Arrays. sz = 1, 2, 3 or 4
    void setUniform(const char* varname, size_t sz, size_t count, float *value);
    void setUniform(unsigned int idVar , size_t sz, size_t count, float *value);
      
    void setUniform(const char* varname, size_t sz, size_t count, int *value);
    void setUniform(unsigned int idVar , size_t sz, size_t count, int *value);

    // dim = 2, 3 or 4
    void setUniformMatrix(const char* varname, size_t dim, size_t count,
                          bool transpose, float* value);
    void setUniformMatrix(unsigned int idVar , size_t dim, size_t count,
                          bool transpose, float* value);
    // Submitting attribute variables

    /** Attribute modificators */
    //@{
    /** Change a scalar vertex attribute */
    void setAttribute(unsigned int index, float v0);
    /** Change a double scalar vertex attribute */
    void setAttribute(unsigned int index, float v0, float v1);
    /** Change a triple scalar vertex attribute */
    void setAttribute(unsigned int index, float v0, float v1, float v2);
    /** Change a quadruple scalar vertex attribute */
    void setAttribute(unsigned int index, float v0, float v1, float v2, float v3);
    /** Change 1 vector vertex attribute 
	dim must be 1, 2, 3 or 4.
     */
    void setAttribute(unsigned int index, unsigned int dim, const float* values);
    //@}
  private:
    void unlink();
    // Tell if this shader program is active or not
    bool _isActive;
    // Tell if shader object is linked
    bool _isLinked;
    // Tell if program must be linked again or not
    bool _isModified;
    // Identificator of vertex/fragment shader programs
    GLhandleARB _programId;
    // List of Vertex/fragment shader subroutines (include unique main subroutine)
    std::list<ShaderObject> _shaders;
  };
}

#endif

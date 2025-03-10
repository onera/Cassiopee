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
#ifndef _CPLOT_SHADER_HPP_
#define _CPLOT_SHADER_HPP_
#include "GL/glew.h"
#include <list>
#include <string>
#include <memory>
#include "VertexShader.h"
#include "TesselationControlShader.hpp"
#include "TesselationEvaluationShader.hpp"
#include "GeomShader.h"
#include "FragmentShader.h"

namespace CPlot
{
  class Shader
  {
  public:
    enum stage
    {
        VERTEX=0, TESSELATION_CONTROL, TESSELATION_EVALUATION, GEOMETRY, FRAGMENT, COMPUTE, END=-1
    };
    std::shared_ptr<VertexShader> vertex_shader;
    std::shared_ptr<TesselationControlShader> tesselation_control_shader;
    std::shared_ptr<TesselationEvaluationShader> tesselation_evaluation_shader;
    std::shared_ptr<GeomShader> geometry_shader;
    std::shared_ptr<FragmentShader> fragment_shader;

    Shader();
    virtual ~Shader();

    void add(const std::shared_ptr<VertexShader>& vshad) 
    {
        _isModified = true; 
        vertex_shader = vshad; 
    }
    void add(const std::shared_ptr<TesselationControlShader>& tcshad) 
    {
        _isModified = true; 
        tesselation_control_shader = tcshad; 
    }
    void add(const std::shared_ptr<TesselationEvaluationShader>& teshad)
    { 
        _isModified = true;
        tesselation_evaluation_shader = teshad; 
    }
    void add(const std::shared_ptr<GeomShader>& gshad)
    {
        _isModified = true;
        geometry_shader = gshad;
    }
    void add(const std::shared_ptr<FragmentShader>& fshad)
    {
        _isModified = true;
        fragment_shader = fshad;
    }
    // void add(const ComputeShader& cshad);


    // Remove a shader from the program (VERTEX, TESSELATION and so.)
    void remove(stage st);
    // Renvoie vrai si contient le proglet passe en parametre...
    bool contain( const ShaderObject& proglet ) const;
    int nb_proglets() const;

    GLhandleARB getProgramId() const { return _programId; }
    bool link();

    bool start(); // Activate the shader. If not linked, link this shader before
    std::string getLinkerLog() const; // Return error log of the linker.
    void end(); // Deactivate the shader. OpenGL calls with the standard pipeline.

    unsigned int getAttributeLocation(const char* attributeName) const;
    unsigned int getUniformLocation  (const char* uniformName  ) const;
    unsigned int getStorageBufferLocation(const char* bufferName  ) const;

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
    void setAttribute(unsigned int index, float v0, float v1,
		      float v2);
    /** Change a quadruple scalar vertex attribute */
    void setAttribute(unsigned int index, float v0, float v1,
		      float v2, float v3);
    /** Change 1 vector vertex attribute 
	dim must be 1, 2, 3 or 4.
     */
    void setAttribute(unsigned int index, unsigned int dim,
		      const float* values);
    //@}

    /** Shader buffer accessors */
    //@{
    /** Initialize a buffer object */
    unsigned int initStorageBuffer( int nfloats, const float* values );
    /** *update data in buffer object */
    void updateStorageBuffer( unsigned int ssbo, int nfloats, const float* values );
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
  };

}

#endif

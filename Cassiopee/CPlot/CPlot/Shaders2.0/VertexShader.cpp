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
# include "GL/glew.h"
# include "ShaderUtil.h"
# include "VertexShader.h"
using namespace CPlot;

//==============================================================================
VertexShader::VertexShader(): ShaderObject()
{
  _shaderId = glCreateShader(GL_VERTEX_SHADER);
  CHECK_GL_ERROR();
}
//==============================================================================
VertexShader::VertexShader(const VertexShader& shader):
  ShaderObject(shader)
{}
//==============================================================================
VertexShader::~VertexShader()
{}
//==============================================================================
const VertexShader&
VertexShader::operator = (const VertexShader& obj)
{
  ShaderObject::operator = (obj);
  return *this;
}

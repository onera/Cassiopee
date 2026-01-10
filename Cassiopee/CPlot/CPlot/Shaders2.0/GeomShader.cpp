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
# include "GeomShader.h"
using namespace CPlot;

//==============================================================================
GeomShader::GeomShader(): ShaderObject()
{
  _shaderId = glCreateShader(GL_GEOMETRY_SHADER);
  CHECK_GL_ERROR();
}
//==============================================================================
GeomShader::GeomShader(const GeomShader& shader):
  ShaderObject(shader)
{}
//==============================================================================
GeomShader::~GeomShader()
{}
//==============================================================================
const GeomShader&
GeomShader::operator = (const GeomShader& obj)
{
  ShaderObject::operator = (obj);
  return *this;
}

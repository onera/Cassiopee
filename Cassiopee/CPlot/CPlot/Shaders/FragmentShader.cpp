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
# include <stdexcept>
# include "GL/glew.h"
# include "FragmentShader.h"
using namespace CPlot;

//==============================================================================
FragmentShader::FragmentShader(): ShaderObject()
{
  _shaderId = glCreateShader(GL_FRAGMENT_SHADER);
  CHECK_GL_ERROR();
}
//==============================================================================
FragmentShader::FragmentShader(const std::string& src):
  ShaderObject()
{
  _shaderId = glCreateShader(GL_FRAGMENT_SHADER);
  bool success = this->compile(src);
  if (not success)
  {
    printf("Warning: CPlot: failed to compile shader.\n");
  }
}
//==============================================================================
FragmentShader::~FragmentShader()
{}
//==============================================================================

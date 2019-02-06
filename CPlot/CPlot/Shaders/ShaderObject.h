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
#ifndef _CPLOT_SHADEROBJECT_HPP_
#define _CPLOT_SHADEROBJECT_HPP_
#include "GL/glew.h"
#include "Def/DefTypes.h"
#include <string>
#include <stdio.h>

namespace CPlot
{
  class ShaderObject
  {
  public:
    ShaderObject();
    ShaderObject(const ShaderObject& obj) = delete;
    virtual ~ShaderObject();

    const ShaderObject& operator = (const ShaderObject& obj) = delete;

    static std::string 
    load(const char* fileName);

    bool compile(const std::string& source);
    std::string getCompilerLog(void) const;
    GLhandleARB getShaderId() const { return _shaderId; }

    unsigned int getAttributeLocation(const char* attributeName) const;
    unsigned int getUniformLocation  (const char* uniformName  ) const;

  protected:
    GLhandleARB _shaderId;

  private:
    void destroy();
    static unsigned E_LONG getFileLength(FILE* ptrFile);    
  };
}

#endif

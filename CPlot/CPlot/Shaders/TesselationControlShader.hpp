/*    
    Copyright 2013-2018 Onera.

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
#ifndef _CPLOT_TESSELATIONCONTROLSHADER_HPP_
#define _CPLOT_TESSELATIONCONTROLSHADER_HPP_
#include "ShaderObject.h"

namespace CPlot
{
  class TesselationControlShader : public ShaderObject
  {
  public:
    TesselationControlShader();
    TesselationControlShader( const std::string& src );
    TesselationControlShader(const TesselationControlShader& shader) = delete;
    virtual ~TesselationControlShader();

    TesselationControlShader& operator = (const TesselationControlShader& obj) = delete;

  };
}

#endif
